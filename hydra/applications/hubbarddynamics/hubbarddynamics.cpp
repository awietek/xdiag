#include <cstdlib>
#include <vector>
#include <utility>
#include <fstream>
#include <iomanip>

#include <lila/all.h>
#include <hydra/all.h>

#include "hubbarddynamics.options.h"

int main(int argc, char* argv[])
{
  using hydra::hilbertspaces::hubbard_qn;
  using hydra::hilbertspaces::Hubbard;
  using hydra::models::HubbardModel;
  using hydra::operators::BondList;
  using hydra::operators::read_bondlist;
  using hydra::dynamics::continued_fraction;

  using namespace lila;

  double t = 1.;
  double U = 4.;
  std::string outfile;
  std::string latticefile;
  int nup = -1;
  int ndown = -1;
  std::string fermiontype;
  int dyniters = 200;
  parse_cmdline(t, U, outfile, latticefile, nup, ndown, fermiontype, dyniters, argc, argv);

  // Open outfile
  std::ofstream of;
  if (outfile != "")
    {
      of.open(outfile, std::ios::out);
      if(of.fail()) 
	{
	  std::cerr << "HubbardThermo Error in opening outfile: " 
		    << "Could not open file with filename ["
		    << outfile << "] given. Abort." << std::endl;
	  exit(EXIT_FAILURE);
	}
    }

  // Parse hoppings from file
  BondList bondlist = read_bondlist(latticefile);
  int n_sites = bondlist.n_sites();
  printf("n_sites: %d\n", n_sites);
  BondList hopping_list = bondlist.bonds_of_type("HUBBARDHOP");
  std::vector<std::pair<int, int>> hoppings;
  for (auto bond : hopping_list)
    hoppings.push_back({bond.sites()[0], bond.sites()[1]});
  auto model = HubbardModel(n_sites, hoppings);

  // Select qns to compute
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};


  printf("Creating Hubbard Hamiltonian for n_upspins=%d, n_downspins=%d...\n", 
	 qn.n_upspins, qn.n_downspins);
  auto hamilton = model.matrix(t, U, qn);
  printf("Done\n");


  printf("Starting ground state eigenvalues Lanczos procedure ...\n");
  double precision = 1e-12;
  int max_iterations = 1000;
  int num_eigenvalue = 0;
  int random_seed = 42;

  auto multiply = [&hamilton](const Vector<double>& v, Vector<double>& w) {
    static int iter=0;
    w = Mult(hamilton, v);
    // printf("iter: %d\n", iter);
    ++iter;
  };
  uint64 dim = hamilton.nrows();
  auto lzs = Lanczos<double, decltype(multiply)>
    (dim, random_seed, max_iterations, precision, num_eigenvalue, multiply);
  Vector<double> eigs = lzs.eigenvalues();
  printf("lzs e %f %f %f %20.18g %20.18g\n", eigs(0), eigs(1), eigs(2), eigs(3), eigs(4));
  printf("Done\n");
   

  // printf("Computing eigenvalues full ED style...\n");
  // auto exeigs = EigenvaluesH(hamilton);
  // printf("ext e %f %f %f %20.18g %20.18g\n", exeigs(0), exeigs(1), exeigs(2), exeigs(3), exeigs(4));
  // printf("Done\n");

  
  printf("Reiterating for ground state...\n");
  auto eigenvectors = lzs.eigenvectors({0});
  Vector<double>& groundstate = eigenvectors[0];
  printf("Done\n");

  LilaPrint(Dot(groundstate, Mult(hamilton, groundstate)));

  printf("Applying creation/annihilation operator...\n");
  Hubbard<uint32> hs(n_sites, qn);
  printf("dim before: %d\n", hs.size());
  Vector<double> dyn_start_state = model.apply_fermion(groundstate, qn, fermiontype, 0);
  double dyn_weight = pow(Norm(dyn_start_state), 2);
  hs = Hubbard<uint32>(n_sites, qn);
  printf("dim after: %d\n", hs.size());  
  printf("Done\n");

  printf("Creating Hubbard Hamiltonian for n_upspins=%d, n_downspins=%d...\n", 
	 qn.n_upspins, qn.n_downspins);
  auto hamilton_dyn = model.matrix(t, U, qn);
  printf("Done\n");


  printf("Starting dynamical Lanczos procedure ...\n");
  precision = -1;
  max_iterations = dyniters;

  auto multiply_dyn = [&hamilton_dyn](const Vector<double>& v, Vector<double>& w) {
    static int iter=0;
    w = Mult(hamilton_dyn, v);
    // printf("dyn iter: %d\n", iter);
    ++iter;
  };

  uint64 dim_dyn = hamilton_dyn.nrows();
  auto lzs_dyn = Lanczos<double, decltype(multiply_dyn)>
    (dim_dyn, random_seed, max_iterations, precision, num_eigenvalue, multiply_dyn);
  lzs_dyn.set_init_state(dyn_start_state);
  
  Vector<double> dyn_eigs = lzs_dyn.eigenvalues();


  // // DEBUG
  // // LilaPrint(lzs_dyn.alphas());
  // // LilaPrint(lzs_dyn.betas());
  // double eta = 0.05;
  // double minz = -5.;
  // double maxz = 5.;
  // int n_z = 10;  
  // for (int k = 0; k <= n_z; ++k)
  //   {
  //     complex thisz = complex(minz + k*(maxz-minz)/n_z, eta);
  //     auto Hinv = thisz * Identity<complex>(hamilton_dyn.nrows()) - Complex(hamilton_dyn);
  //     Invert(Hinv);
  //     printf("z: %f + 1j*%f\n", std::real(thisz), std::imag(thisz));
  //     auto state_c = Complex(dyn_start_state);
  //     auto tmp = Mult(Hinv, state_c);
  //     LilaPrint(Dot(state_c, tmp));
  //     LilaPrint(continued_fraction(thisz, lzs_dyn.alphas(), lzs_dyn.betas()) * dyn_weight);
  //   }
    
  // Write to outfile
  if (outfile != "")
    {
      std::stringstream line;
      line << std::setprecision(20);
      line << "# gs energy: " << eigs(0) << "\n";
      line << "# weight: " << dyn_weight << "\n";
      line << "# alphas betas\n";
      of << line.str();
      line.str("");

      auto alphas = lzs_dyn.alphas();
      auto betas = lzs_dyn.betas();
      assert(alphas.size() == betas.size());
      
      for (int idx = 0; idx < alphas.size(); ++idx)
	{
	  line << alphas(idx) << " " << betas(idx) << "\n";
	  // std::cout << line.str();
	  of << line.str();
	  line.str("");
	}
    }
  return EXIT_SUCCESS;
}
