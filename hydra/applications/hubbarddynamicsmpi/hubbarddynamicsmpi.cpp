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
#include "hubbarddynamicaliterationsmpi.h"

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
  using hydra::hilbertspaces::hubbard_qn;
  using hydra::hilbertspaces::Hubbard;
  using hydra::models::HubbardModelMPI;
  using namespace hydra::operators;
  using hydra::dynamics::continued_fraction;
  using namespace lila;

  MPI_Init(&argc, &argv); 
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  if (mpi_rank == 0) printf("Using %d MPI tasks\n", mpi_size);

  // Get input parameters
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
  int verbosity = 1;
  double deflationtol = 1e-8;

  parse_cmdline(outfile, latticefile, couplingfile, corrfile, nup, ndown, 
		fermiontype, algorithm, precision, iters, dynprecision, 
		dyniters, verbosity, deflationtol, argc, argv);

  // Check if valid algorithm is defined
  if (algorithm == "") algorithm = "lanczos";
  if (!((algorithm == "lanczos") || (algorithm == "bandlanczos")))
    {
      if (mpi_rank == 0) 
	printf("Unknown algorithm: %s (expected lanczos or bandlanczos)\n", 
	       algorithm.c_str());
      MPI_Finalize();
      return EXIT_FAILURE;
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
	  if (mpi_rank == 0)
	    {
	      std::cerr << "Error in opening outfile: " 
			<< "Could not open file with filename ["
			<< outfile << "] given. Abort." << std::endl;
	      MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	}
    }


  // Parse bondlist and couplings from file
  BondList bondlist = read_bondlist(latticefile);
  BondList hopping_list = bondlist.bonds_of_type("HUBBARDHOP");
  BondList interaction_list = bondlist.bonds_of_type("HUBBARDV");
  if (couplingfile == "") couplingfile = latticefile;
  Couplings couplings = read_couplings(couplingfile);

  if ((verbosity >= 1) && (mpi_rank == 0))
    {
      for (auto bond : hopping_list)
	printf("hopping %s %d %d\n", bond.coupling().c_str(), 
	       bond.sites()[0], bond.sites()[1]);
      for (auto bond : interaction_list)
	printf("interaction %s %d %d\n", bond.coupling().c_str(), 
	       bond.sites()[0], bond.sites()[1]);
      for (auto c : couplings)
	printf("coupling %s %f %fj\n", c.first.c_str(), std::real(c.second), std::imag(c.second));
    }

  int n_sites = bondlist.n_sites();
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};


  // Parse the corrfile
  std::vector<std::pair<int, int>> correlation_list;
  if (corrfile == "") correlation_list.push_back({0, 0});
  else
    {
      std::ifstream ifs(corrfile);
      if(ifs.fail()) 
	{
	  if (mpi_rank == 0)
	    std::cerr << "Error in opening corrfile: " 
		      << "Could not open file with filename ["
		      << outfile << "] given. Abort." << std::endl;
	  MPI_Abort(MPI_COMM_WORLD, 2);
	}
      int a, b;
      while (ifs >> a >> b) correlation_list.push_back({a, b});
    }

  for (auto corr : correlation_list)
    {
      if ((corr.first >= n_sites) || (corr.second >= n_sites))
	{
	  if (mpi_rank == 0) printf("Invalid correlation: %d %d\n", corr.first, corr.second);
	  MPI_Abort(MPI_COMM_WORLD, 3);
	}
      if (verbosity >= 2)  
	if (mpi_rank == 0) printf("corr %d %d\n", corr.first, corr.second);
    }


  // Create infrastructure for Hubbard model
  if (mpi_rank == 0) 
    printf("Creating Hubbard model for n_upspins=%d, n_downspins=%d...\n",
	   qn.n_upspins, qn.n_downspins);
  double t1 = MPI_Wtime();
  auto model = HubbardModelMPI<uint32>(bondlist, couplings, qn);
  double t2 = MPI_Wtime();
  if (mpi_rank == 0) printf("time init: %3.4f\n", t2-t1); 
  


  // Get the ground state
  if (mpi_rank == 0) 
    printf("Starting ground state eigenvalues Lanczos procedure ...\n");
  int num_eigenvalue = 0;
  int random_seed = 42 + 1234567*mpi_rank;

  auto multiply = [&model, &mpi_rank, &verbosity](const VectorMPI<double>& v, 
				      VectorMPI<double>& w) {
    static int iter=0;
    double t1 = MPI_Wtime();
    bool verbose = ((iter == 0) && verbosity >=1);
    model.apply_hamiltonian(v, w, verbose);
    double t2 = MPI_Wtime();
    if ((verbosity >= 2) && (mpi_rank == 0))
      printf("iter: %d, time MVM: %3.4f\n", iter, t2-t1); 
    ++iter;
  };
  MPI_Barrier(MPI_COMM_WORLD);
  uint64 dim = model.dim();
  if (mpi_rank == 0) printf("dim %s\n", FormatWithCommas(dim).c_str()); 
  MPI_Barrier(MPI_COMM_WORLD);
  uint64 local_dim = model.local_dim();
  auto lzs = Lanczos<double, decltype(multiply), VectorMPI<double>>
    (local_dim, random_seed, iters, precision, num_eigenvalue, multiply);
  Vector<double> eigs = lzs.eigenvalues();

  if (mpi_rank == 0) printf("Done\n");
  

  if (mpi_rank == 0) printf("Reiterating for ground state...\n");
  auto eigenvectors = lzs.eigenvectors({0});
  VectorMPI<double>& groundstate = eigenvectors[0];
  if (mpi_rank == 0) printf("Done\n");

  // Write gs energy to outfile
  if (outfile != "")
    {
      if (mpi_rank == 0)
	{
	  std::stringstream line;
	  line << std::setprecision(20);
	  line << "# gs energy: " << eigs(0) << "\n";
	  of << line.str();
	  line.str("");
	}
    }

  if (algorithm == "lanczos")
    {
      for (auto corr : correlation_list)
	for (auto ftype : ftype_list)
	  {  
	    int s1 = corr.first;
	    int s2 = corr.second;
	    auto res = hydra::hubbard_dynamical_iterations_lanczos_mpi
	      (model, groundstate, s1, s2, ftype, dyniters, dynprecision, verbosity);

	    // Write to outfile
	    if (mpi_rank == 0)
	      {
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
	  auto res = hydra::hubbard_dynamical_iterations_bandlanczos_mpi
	    (model, groundstate, sites, ftype, dyniters, dynprecision, 
	     verbosity, deflationtol);
	  auto tmat = res.tmatrix;
	  auto overlaps = res.overlaps;
	  std::stringstream line;
	  if (mpi_rank == 0)
	    {
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
	}
    }  // algorithm bandlanczos

  MPI_Finalize();
  return EXIT_SUCCESS;
}
