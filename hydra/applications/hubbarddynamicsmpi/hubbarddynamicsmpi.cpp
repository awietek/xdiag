#include <cstdlib>
#include <vector>
#include <utility>
#include <fstream>
#include <iomanip>
#include <unistd.h>

#include <mpi.h>

#include <lila/allmpi.h>
#include <hydra/allmpi.h>

#include "hubbarddynamicsmpi.options.h"

#include <iomanip>
#include <locale>

template<class T>
std::string FormatWithCommas(T value)
{
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << value;
    return ss.str();
}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv); 
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  using hydra::hilbertspaces::hubbard_qn;
  using hydra::hilbertspaces::Hubbard;
  using hydra::models::HubbardModelMPI;
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
	  if (mpi_rank == 0)
	    {
	      std::cerr << "HubbardDynamicsMPI Error in opening outfile: " 
			<< "Could not open file with filename ["
			<< outfile << "] given. Abort." << std::endl;
	    }
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}
    }

  // Parse hoppings from file
  BondList bondlist = read_bondlist(latticefile);
  int n_sites = bondlist.n_sites();
  BondList hopping_list = bondlist.bonds_of_type("HUBBARDHOP");
  std::vector<std::pair<int, int>> hoppings;
  for (auto bond : hopping_list)
    hoppings.push_back({bond.sites()[0], bond.sites()[1]});

  // select qns to compute
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};

  double t1 = MPI_Wtime();
  auto model = HubbardModelMPI<uint32>(n_sites, hoppings, qn);
  double t2 = MPI_Wtime();
  if (mpi_rank == 0) printf("time init: %3.4f\n", t2-t1); 
  

  // if (mpi_rank == 0) printf("Sleeping after creating model ...\n");
  // sleep(20);


  if (mpi_rank == 0) printf("Starting ground state eigenvalues Lanczos procedure ...\n");
  double precision = 1e-12;
  int max_iterations = 1000;
  int num_eigenvalue = 0;
  int random_seed = 42 + 1234567*mpi_rank;

  auto multiply = [&model, &t, &U, &mpi_rank](const VectorMPI<double>& v, VectorMPI<double>& w) {
    static int iter=0;
    double t1 = MPI_Wtime();
    model.apply_hamiltonian(t, U, v, w);
    double t2 = MPI_Wtime();
    if (mpi_rank == 0) printf("time MVM: %3.4f\n", t2-t1); 
    // printf("iter: %d\n", iter);
    ++iter;
  };
  MPI_Barrier(MPI_COMM_WORLD);
  uint64 dim = model.dim();
  if (mpi_rank == 0) printf("dim %s\n", FormatWithCommas(dim).c_str()); 
  MPI_Barrier(MPI_COMM_WORLD);
  uint64 local_dim = model.local_dim();
  printf("[%d] ld: %d\n", mpi_rank, local_dim);
  // uint64 dim = model.dim();
  // printf("[%d] d: %d\n", mpi_rank, dim);
  auto lzs = Lanczos<double, decltype(multiply), VectorMPI<double>>
    (local_dim, random_seed, max_iterations, precision, num_eigenvalue, multiply);
  Vector<double> eigs = lzs.eigenvalues();
  if (mpi_rank == 0) LilaPrint(eigs);
  if (mpi_rank == 0) printf("Done\n");
  

  if (mpi_rank == 0) printf("Reiterating for ground state...\n");
  auto eigenvectors = lzs.eigenvectors({0});
  VectorMPI<double>& groundstate = eigenvectors[0];
  if (mpi_rank == 0) printf("Done\n");

  auto w = groundstate;
  multiply(groundstate, w);
  auto e0 = Dot(groundstate, w);
  if (mpi_rank == 0) LilaPrint(e0);

  // LilaPrint(Dot(groundstate, Mult(hamilton, groundstate)));

  // printf("Applying creation/annihilation operator...\n");
  // Hubbard<uint32> hs(n_sites, qn);
  // printf("dim before: %d\n", hs.size());
  // Vector<double> dyn_start_state = model.apply_fermion(groundstate, qn, fermiontype, 0);
  // double dyn_weight = pow(Norm(dyn_start_state), 2);
  // hs = Hubbard<uint32>(n_sites, qn);
  // printf("dim after: %d\n", hs.size());  
  // printf("Done\n");

  // printf("Creating Hubbard Hamiltonian for n_upspins=%d, n_downspins=%d...\n", 
  // 	 qn.n_upspins, qn.n_downspins);
  // auto hamilton_dyn = model.matrix(t, U, qn);
  // printf("Done\n");


  // printf("Starting dynamical Lanczos procedure ...\n");
  // precision = -1;
  // max_iterations = dyniters;

  // auto multiply_dyn = [&hamilton_dyn](const Vector<double>& v, Vector<double>& w) {
  //   static int iter=0;
  //   w = Mult(hamilton_dyn, v);
  //   // printf("dyn iter: %d\n", iter);
  //   ++iter;
  // };

  // uint64 dim_dyn = hamilton_dyn.nrows();
  // auto lzs_dyn = Lanczos<double, decltype(multiply_dyn)>
  //   (dim_dyn, random_seed, max_iterations, precision, num_eigenvalue, multiply_dyn);
  // lzs_dyn.set_init_state(dyn_start_state);
  
  // Vector<double> dyn_eigs = lzs_dyn.eigenvalues();


  // // // DEBUG
  // // // LilaPrint(lzs_dyn.alphas());
  // // // LilaPrint(lzs_dyn.betas());
  // // double eta = 0.05;
  // // double minz = -5.;
  // // double maxz = 5.;
  // // int n_z = 10;  
  // // for (int k = 0; k <= n_z; ++k)
  // //   {
  // //     complex thisz = complex(minz + k*(maxz-minz)/n_z, eta);
  // //     auto Hinv = thisz * Identity<complex>(hamilton_dyn.nrows()) - Complex(hamilton_dyn);
  // //     Invert(Hinv);
  // //     printf("z: %f + 1j*%f\n", std::real(thisz), std::imag(thisz));
  // //     auto state_c = Complex(dyn_start_state);
  // //     auto tmp = Mult(Hinv, state_c);
  // //     LilaPrint(Dot(state_c, tmp));
  // //     LilaPrint(continued_fraction(thisz, lzs_dyn.alphas(), lzs_dyn.betas()) * dyn_weight);
  // //   }
    
  // // Write to outfile
  // if (outfile != "")
  //   {
  //     std::stringstream line;
  //     line << std::setprecision(20);
  //     line << "# gs energy: " << eigs(0) << "\n";
  //     line << "# weight: " << dyn_weight << "\n";
  //     line << "# alphas betas\n";
  //     of << line.str();
  //     line.str("");

  //     auto alphas = lzs_dyn.alphas();
  //     auto betas = lzs_dyn.betas();
  //     assert(alphas.size() == betas.size());
      
  //     for (int idx = 0; idx < alphas.size(); ++idx)
  // 	{
  // 	  line << alphas(idx) << " " << betas(idx) << "\n";
  // 	  // std::cout << line.str();
  // 	  of << line.str();
  // 	  line.str("");
  // 	}
  //   }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
