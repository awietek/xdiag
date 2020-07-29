#include <cstdlib>

#include <mpi.h>

#include <lila/allmpi.h>
#include <hydra/allmpi.h>
#include <lime/all.h>

#include <stdio.h>

lila::LoggerMPI lg;

#include "tjdynamicsmpi.options.h"
#include "tjdynamicaliterationsmpi.h"

template <class coeff_t>
void run_real_complex(std::string real_complex,
    hydra::all::BondList bondlist,
    hydra::all::Couplings couplings,
    hydra::all::BondList corr_bondlist,
    hydra::all::hubbard_qn qn,
    double precision, 
    int iters, std::string algorithm,
    double dynprecision, int dyniters, 
    double deflationtol, int lobpcgbands,
    int verbosity, int seed,
    std::string outfile)
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  FileH5 file;
  if (mpi_rank == 0) file = lime::FileH5(outfile, "w!");

  lg.out(1, "Creating {} t-J model for n_upspins={}, n_downspins={}...\n",
      real_complex, qn.n_upspins, qn.n_downspins);
  lg.out(1, "Using {} MPI tasks\n", mpi_size);
  
  // Create tJ Hamiltonian 

  double t1 = MPI_Wtime();
  auto H = TJModelMPI<coeff_t>(bondlist, couplings, qn);
  double t2 = MPI_Wtime();
  lg.out(1, "time init: {} secs\n", t2-t1); 


  // Define Hamilton multiplication function
  int iter = 0;
  auto multiply_H = 
    [&H, &iter, verbosity](const VectorMPI<coeff_t>& v, VectorMPI<coeff_t>& w) 
    {
      bool verbose = (iter==0) && (verbosity > 0);
      double t1 = MPI_Wtime();
      H.apply_hamiltonian(v, w, verbose);
      double t2 = MPI_Wtime();
      lg.out(2, "iter: {}, time MVM: {}\n", iter, t2-t1); 
      ++iter;
    };


  ///////////////////////////////
  // Get the ground state
  auto local_dim = H.local_dim();
  VectorMPI<coeff_t> groundstate(local_dim);
  Vector<double> eigs;

  int random_seed = seed + 1234567*mpi_rank;
  normal_dist_t<coeff_t> dist(0., 1.);
  normal_gen_t<coeff_t> gen(dist, random_seed);

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
      std::vector<VectorMPI<coeff_t>> vs;
      for (int i=0; i<lobpcgbands; ++i)
	{
	  VectorMPI<coeff_t> v(local_dim);
	  Random(v, gen);
	  vs.push_back(v);
	}

      // Run LOBPCG
      auto res = Lobpcg(multiply_H, vs, precision, iters);
      eigs = res.eigenvalues;
      groundstate = res.eigenvectors[0];
    }

  coeff_t norm = Dot(groundstate, groundstate);
  groundstate /= sqrt(norm);
  lg.out(1, "Done\n");
  lg.out(1, "Ground state energy: {}\n", eigs(0));
  if (mpi_rank == 0)
  {
    file["GroundStateEnergy"] = eigs(0);
  }


  ///////////////////////////////
  // Get dynamics
  if (algorithm == "lanczos")
    {
      auto corrs_sz = corr_bondlist.bonds_of_type("SZ");
      for (auto corr : corrs_sz)
	{
	  // Get the dynamical weight and Lanczos T-Matrix
	  assert(corr.sites().size() == 2);
	  int s1 = corr.sites()[0];
	  int s2 = corr.sites()[1];
	  auto res = hydra::tj_dynamical_iterations_lanczos_mpi<coeff_t>
	    (H, groundstate, s1, s2, "sz", dyniters, dynprecision);
	  
	  // Dump to outfile
	  std::stringstream ss;
	  ss << "sz_" << s1 << "_" << s2;
      std::string label = ss.str();
      if (mpi_rank == 0)
      {
      file[label + std::string("_Weight")] = res.dyn_weight;
      file[label + std::string("_Alphas")] = res.alphas;
      file[label + std::string("_Betas")] = res.betas;
      file[label + std::string("_Eigenvalues")] = res.eigenvalues;
  }
	} 
    }
  if (mpi_rank == 0) {
	file.close();
  }
 
}

int main(int argc, char* argv[])
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;

  MPI_Init(&argc, &argv); 

  // Get input parameters
  std::string outfile;
  std::string latticefile;
  std::string couplingfile;
  std::string corrfile;
  int nup = -1;
  int ndown = -1;
  std::string algorithm;
  double precision = 1e-12;
  int iters = 1000;
  double dynprecision = 1e-12;
  int dyniters = 1000;
  int seed = 1;
  int verbosity = 1;
  double deflationtol = 1e-8;
  int lobpcgbands = 1;
  parse_cmdline(outfile, latticefile, couplingfile, corrfile, nup, ndown, 
		 algorithm, precision, iters, dynprecision, 
		dyniters, verbosity, deflationtol, lobpcgbands, seed, argc, argv);
  lg.set_verbosity(verbosity);  
  check_if_files_exists({latticefile, couplingfile, corrfile});
  check_if_contained_in(algorithm, {"lanczos", "bandlanczos"});


  // Parse bondlist, couplings, and corr_bondlist from file
  BondList bondlist = read_bondlist(latticefile);
  Couplings couplings = read_couplings(couplingfile);
  BondList corr_bondlist = read_bondlist(corrfile);


  int n_sites = bondlist.n_sites();
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};


  if (couplings.all_real())
    run_real_complex<double>(std::string("REAL"), bondlist,
        couplings, corr_bondlist, qn, precision, iters, algorithm,
        dynprecision, dyniters, deflationtol, lobpcgbands,
        verbosity, seed, outfile);
  else
    run_real_complex<lila::complex>(std::string("COMPLEX"), bondlist,
        couplings, corr_bondlist, qn, precision, iters, algorithm,
        dynprecision, dyniters, deflationtol, lobpcgbands,
        verbosity, seed, outfile);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
  
