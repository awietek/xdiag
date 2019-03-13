#include <cstdlib>
#include <vector>
#include <utility>
#include <fstream>
#include <iomanip>

#include <lila/all.h>
#include <hydra/all.h>

#include "hubbarddynamics.options.h"
#include "hubbarddynamicaliterations.h"

int main(int argc, char* argv[])
{
  using hydra::hilbertspaces::hubbard_qn;
  using hydra::hilbertspaces::Hubbard;
  using hydra::models::HubbardModel;
  using namespace hydra::operators;
  using Clock = std::chrono::high_resolution_clock;
  using secs = std::chrono::duration<float>;
  using namespace lila;

  std::string outfile;
  std::string latticefile;
  std::string couplingfile;
  std::string corrfile;
  int nup = -1;
  int ndown = -1;
  std::string fermiontype;
  std::string algorithm;
  double precision = 1e-12;
  int iters = 1000;
  double dynprecision = 1e-12;
  int dyniters = 1000;
  parse_cmdline(outfile, latticefile, couplingfile, corrfile, nup, ndown, 
		fermiontype, algorithm, precision, iters, dynprecision, 
		dyniters, argc, argv);

  // Check if valid algorithm is defined
  if (algorithm == "") algorithm = "lanczos";
  if (!((algorithm == "lanczos") || (algorithm == "bandlanczos")))
    {
      printf("Unknown algorithm: %s (expected lanczos or bandlanczos)\n", 
	     algorithm.c_str());
      exit(EXIT_FAILURE);
    } 

  // Set fermiontype list
  std::vector<std::string> ftype_list;
  if (fermiontype != "") ftype_list.push_back(fermiontype);
  else
    {
      ftype_list.push_back("cdagdn");
      ftype_list.push_back("cdn");
    }

  // Open outfile
  std::ofstream of;
  if (outfile != "")
    {
      of.open(outfile, std::ios::out);
      if(of.fail()) 
	{
	  std::cerr << "Error in opening outfile: " 
		    << "Could not open file with filename ["
		    << outfile << "] given. Abort." << std::endl;
	  exit(EXIT_FAILURE);
	}
    }

  // Parse bondlist and couplings from file
  BondList bondlist = read_bondlist(latticefile);
  BondList hopping_list = bondlist.bonds_of_type("HUBBARDHOP");
  for (auto bond : hopping_list)
    printf("hopping %s %d %d\n", bond.coupling().c_str(), 
	   bond.sites()[0], bond.sites()[1]);
  BondList interaction_list = bondlist.bonds_of_type("HUBBARDV");
  for (auto bond : interaction_list)
    printf("interaction %s %d %d\n", bond.coupling().c_str(), 
	   bond.sites()[0], bond.sites()[1]);
  if (couplingfile == "") couplingfile = latticefile;
  Couplings couplings = read_couplings(couplingfile);
  for (auto c : couplings)
    printf("coupling %s %f %fj\n", c.first.c_str(), std::real(c.second), std::imag(c.second));
  int n_sites = bondlist.n_sites();
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};
  auto model = HubbardModel(bondlist, couplings, qn);

  // Parse the corrfile
  std::vector<std::pair<int, int>> correlation_list;
  if (corrfile == "") correlation_list.push_back({0, 0});
  else
    {
      std::ifstream ifs(corrfile);
      if(ifs.fail()) 
	{
	  std::cerr << "Error in opening corrfile: " 
		    << "Could not open file with filename ["
		    << outfile << "] given. Abort." << std::endl;
	  exit(EXIT_FAILURE);
	}
      int a, b;
      while (ifs >> a >> b) correlation_list.push_back({a, b});
    }
  for (auto corr : correlation_list)
    {
      if ((corr.first >= n_sites) || (corr.second >= n_sites))
	{
	  printf("Invalid correlation: %d %d\n", corr.first, corr.second);
	  exit(EXIT_FAILURE);
	}
      printf("measure corr %d %d\n", corr.first, corr.second);
    }


  printf("Starting ground state eigenvalues Lanczos procedure ...\n");

  int num_eigenvalue = 0;
  int random_seed = 42;

  auto multiply = [&model](const Vector<double>& v, Vector<double>& w) {
    static int iter=0;
    auto t1 = Clock::now();
    model.apply_hamiltonian(v, w);
    auto t2 = Clock::now();
    printf("iter: %d, time MVM: %3.4f\n", iter, secs(t2-t1).count()); 
    ++iter;
  };
  
  int64 dim = model.dim();
  auto lzs = Lanczos<double, decltype(multiply)>
    (dim, random_seed, iters, precision, num_eigenvalue, multiply);
  Vector<double> eigs = lzs.eigenvalues();
  printf("lzs e %20.18g\n", eigs(0));
  printf("Done\n");
  
  printf("Reiterating for ground state...\n");
  auto eigenvectors = lzs.eigenvectors({0});
  Vector<double>& groundstate = eigenvectors[0];
  printf("Done\n");

  // Write gs energy to outfile
  if (outfile != "")
    {
      std::stringstream line;
      line << std::setprecision(20);
      line << "# gs energy: " << eigs(0) << "\n";
      of << line.str();
      line.str("");
    }


  if (algorithm == "lanczos")
    {
      for (auto corr : correlation_list)
	for (auto ftype : ftype_list)
	  {  
	    int s1 = corr.first;
	    int s2 = corr.second;
	    auto res = hydra::hubbard_dynamical_iterations_lanczos
	      (model, groundstate, s1, s2, ftype, dyniters, dynprecision);

	    // Write to outfile
	    std::stringstream line;
	    line << "# ftype: " << ftype << ", s1: " << s1 << ", s2: " 
		 << s2 << "\n";
	    line << "# weight: " << std::setprecision(20) 
		 << res.dyn_weight << "\n";
	    line << "# alphas betas\n";
	    of << line.str();
	    line.str("");

	    auto alphas = res.alphas;
	    auto betas = res.betas;
	    assert(alphas.size() == betas.size());
	    for (int idx = 0; idx < alphas.size(); ++idx)
	      {
		line << std::setprecision(20) 
		     << alphas(idx) << " " << betas(idx) << "\n";
		of << line.str();
		line.str("");
	      }	      
	  }
    }  // algorithm lanczos
  else if (algorithm == "bandlanczos")
    {
      std::vector<int> sites;
      for (auto corr : correlation_list)
	{
	  int s1 = corr.first;
	  int s2 = corr.second;
	  if (std::find(sites.begin(), sites.end(), s1) == sites.end())
	    sites.push_back(s1);
	  if (std::find(sites.begin(), sites.end(), s2) == sites.end())
	    sites.push_back(s2);
	}
      for (auto ftype : ftype_list)
	{  
	  auto res = hydra::hubbard_dynamical_iterations_bandlanczos
	    (model, groundstate, sites, ftype, dyniters, dynprecision);
	  auto tmat = res.tmatrix;
	  auto overlaps = res.overlaps;
	  std::stringstream line;
	  line << "# ftype: " << ftype << "\n";
	  line << "# sites: ";
	  for (auto site : sites)
	    line << site << " ";
	  line << "\n";
	  line << "# tmatrix\n";
	  of << line.str();

	  // Write the tmatrix
	  line.str("");
	  for (auto i : tmat.rows())
	    {
	      for (auto j : tmat.cols())
		line << std::setprecision(20) << tmat(i, j) << " ";
	      line << "\n";
	    }
	  of << line.str();

	  // Write the overlaps
	  line.str("");
	  line << "# overlaps\n";
	  of << line.str();
	  line.str("");
	  for (auto i : overlaps.rows())
	    {
	      for (auto j : overlaps.cols())
		line << std::setprecision(20) << overlaps(i, j) << " ";
	      line << "\n";
	    }
	  of << line.str();

	}
    }  // algorithm bandlanczos
 
  return EXIT_SUCCESS;
}
