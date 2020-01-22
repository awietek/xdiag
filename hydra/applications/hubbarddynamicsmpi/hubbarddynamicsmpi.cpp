#include <cstdlib>

#include <mpi.h>

#include <lila/allmpi.h>
#include <hydra/allmpi.h>
#include <lime/all.h>

lila::LoggerMPI lg;

#include "hubbarddynamicsmpi.options.h"
#include "hubbarddynamicaliterationsmpi.h"

int main(int argc, char* argv[])
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;

  MPI_Init(&argc, &argv); 
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);


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
  int lobpcgbands = 1;
  parse_cmdline(outfile, latticefile, couplingfile, corrfile, nup, ndown, 
		fermiontype, algorithm, precision, iters, dynprecision, 
		dyniters, verbosity, deflationtol, lobpcgbands, argc, argv);
  lg.set_verbosity(verbosity);  
  check_if_files_exists({latticefile, couplingfile, corrfile});
  check_if_contained_in(algorithm, {"lanczos", "bandlanczos"});
  auto dumper = lime::MeasurementsH5((mpi_rank == 0) ? outfile : "");
  lg.out(1, "Using {} MPI tasks\n", mpi_size);


  // Parse bondlist, couplings, and corr_bondlist from file
  BondList bondlist = read_bondlist(latticefile);
  Couplings couplings = read_couplings(couplingfile);
  BondList corr_bondlist = read_bondlist(corrfile);


  int n_sites = bondlist.n_sites();
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};


  // Create Hubbard Hamiltonian
  lg.out(1, "Creating Hubbard model for n_upspins={}, n_downspins={}...\n",
	     qn.n_upspins, qn.n_downspins);
  double t1 = MPI_Wtime();
  auto H = HubbardModelMPI<double,uint32>(bondlist, couplings, qn);
  double t2 = MPI_Wtime();
  lg.out(1, "time init: {} secs\n", t2-t1); 


  // Define Hamilton multiplication function
  int iter = 0;
  auto multiply_H = 
    [&H, &iter](const VectorMPI<double>& v, VectorMPI<double>& w) 
    {
      double t1 = MPI_Wtime();
      H.apply_hamiltonian(v, w, false);
      double t2 = MPI_Wtime();
      lg.out(2, "iter: {}, time MVM: {}\n", iter, t2-t1); 
      ++iter;
    };


  ///////////////////////////////
  // Get the ground state
  uint64 local_dim = H.local_dim();
  VectorMPI<double> groundstate(local_dim);
  Vector<double> eigs;

  int random_seed = 42 + 1234567*mpi_rank;
  normal_dist_t<double> dist(0., 1.);
  normal_gen_t<double> gen(dist, random_seed);

  // Lanczos algorithm for ground state (unstable)
  if (lobpcgbands == 0)
    {
      lg.out(1, "Starting ground state Lanczos...\n");
      auto res = LanczosEigenvectors(multiply_H, groundstate, gen, true, precision, {0}, "Ritz");
      eigs = res.eigenvalues;
      groundstate = res.vectors[0];
    }

  // LOBPCG algorithm for ground state (prefered)
  else
    {
      lg.out(1, "Starting ground state LOBPCG algorithm...\n");
      lg.out(1, "Using {} bands\n", lobpcgbands);
      
      // Create random start vectors
      std::vector<VectorMPI<double>> vs;
      for (int i=0; i<lobpcgbands; ++i)
	{
	  VectorMPI<double> v(local_dim);
	  Random(v, gen);
	  vs.push_back(v);
	}

      // Run LOBPCG
      auto res = Lobpcg(multiply_H, vs, precision, iters);
      eigs = res.eigenvalues;
      groundstate = res.eigenvectors[0];
    }

  double norm = Dot(groundstate, groundstate);
  groundstate /= sqrt(norm);
  lg.out(1, "Done\n");
  lg.out(1, "Ground state energy: {}\n", eigs(0));
  lg.out(1, "Ground state norm  : {}\n", norm);
  dumper["GroundStateEnergy"] << eigs(0);
  dumper.dump();

  ///////////////////////////////
  // Get dynamics
  if (algorithm == "lanczos")
    {
      auto corrs_cdagdncdn = corr_bondlist.bonds_of_type("CDAGDNCDN");
      for (auto corr : corrs_cdagdncdn)
	{
	  // Get the dynamical weight and Lanczos T-Matrix
	  assert(corr.sites().size() == 2);
	  int s1 = corr.sites()[0];
	  int s2 = corr.sites()[1];
	  auto res = hydra::hubbard_dynamical_iterations_lanczos_mpi
	    (H, groundstate, s1, s2, "cdagdn", dyniters, dynprecision);
	  
	  // Dump to outfile
	  std::stringstream ss;
	  ss << "CDAGDNCDN_" << s1 << "_" << s2;
	  std::string label = ss.str();
	  dumper[label + std::string("_Weight")] << res.dyn_weight;
	  dumper[label + std::string("_Alphas")] << res.alphas;
	  dumper[label + std::string("_Betas")] << res.betas;
	  dumper[label + std::string("_Eigenvalues")] << res.eigenvalues;
	  dumper.dump();
	} 
    }

  // else if (algorithm == "bandlanczos")
  //   {
  //     std::vector<int> sites;
  //     for (auto corr : correlation_list)
  // 	{
  // 	  int s1 = corr.first;
  // 	  int s2 = corr.second;
  // 	  if (std::find(sites.begin(), sites.end(), s1) == sites.end())
  // 	    sites.push_back(s1);
  // 	  if (std::find(sites.begin(), sites.end(), s2) == sites.end())
  // 	    sites.push_back(s2);
  // 	}
  //     for (auto ftype : ftype_list)
  // 	{  
  // 	  auto res = hydra::hubbard_dynamical_iterations_bandlanczos_mpi
  // 	    (model, groundstate, sites, ftype, dyniters, dynprecision, 
  // 	     verbosity, deflationtol);
  // 	  auto tmat = res.tmatrix;
  // 	  auto overlaps = res.overlaps;
  // 	  std::stringstream line;
  // 	  if (mpi_rank == 0)
  // 	    {
  // 	      line << "# ftype: " << ftype << "\n";
  // 	      line << "# sites: ";
  // 	      for (auto site : sites)
  // 		line << site << " ";
  // 	      line << "\n";
  // 	      line << "# tmatrix\n";
  // 	      of << line.str();

  // 	      // Write the tmatrix
  // 	      line.str("");
  // 	      for (auto i : tmat.rows())
  // 		{
  // 		  for (auto j : tmat.cols())
  // 		    line << std::setprecision(20) << tmat(i, j) << " ";
  // 		  line << "\n";
  // 		}
  // 	      of << line.str();

  // 	      // Write the overlaps
  // 	      line.str("");
  // 	      line << "# overlaps\n";
  // 	      of << line.str();
  // 	      line.str("");
  // 	      for (auto i : overlaps.rows())
  // 		{
  // 		  for (auto j : overlaps.cols())
  // 		    line << std::setprecision(20) << overlaps(i, j) << " ";
  // 		  line << "\n";
  // 		}
  // 	      of << line.str();
  // 	    }
  // 	}
  //   }  // algorithm bandlanczos
  
  MPI_Finalize();
  return EXIT_SUCCESS;
}
