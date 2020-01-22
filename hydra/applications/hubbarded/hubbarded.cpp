#include <cstdlib>

#include <mpi.h>

#include <lila/allmpi.h>
#include <hydra/allmpi.h>
#include <lime/all.h>

lila::LoggerMPI lg;

#include "hubbarded.options.h"

int main(int argc, char* argv[])
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;

  MPI_Init(&argc, &argv); 
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  std::string outfile;
  std::string latticefile;
  std::string couplingfile;
  int nup = -1;
  int ndown = -1;
  double precision = 1e-12;
  int neigenvalue = 0;
  int iters = 1000;
  int verbosity = 1;
  int seed = 1;
  parse_cmdline(outfile, latticefile, couplingfile, nup, ndown, precision, 
		neigenvalue, iters, verbosity, seed, argc, argv);

  lg.set_verbosity(verbosity);  
  
  auto dumper = lime::MeasurementsH5((mpi_rank == 0) ? outfile : "");
  check_if_files_exists({latticefile, couplingfile});

  // Create Hamiltonian
  BondList bondlist = read_bondlist(latticefile);
  Couplings couplings = read_couplings(couplingfile);

  // Create infrastructure for Hubbard model
  int n_sites = bondlist.n_sites();
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};
      
  lg.out(1, "Creating Hubbard model for n_upspins={}, n_downspins={}...\n",
	     qn.n_upspins, qn.n_downspins);
  auto H = HubbardModelMPI<double>(bondlist, couplings, qn);

  // Define multiplication function
  int iter = 0;
  auto multiply_H = 
    [&H, &iter](const VectorMPI<double>& v, VectorMPI<double>& w) 
    {
      H.apply_hamiltonian(v, w);
      lg.out(2, "iter: {}\n", iter); 
      ++iter;
    };

  // Create normal distributed random start state
  VectorMPI<double> startstate(H.dim());
  normal_dist_t<double> dist(0., 1.);
  normal_gen_t<double> gen(dist, seed);
  Random(startstate, gen, true);
  Normalize(startstate);  

  // Run Lanczos
  auto res = LanczosEigenvalues(multiply_H, startstate, precision,
				neigenvalue, "Eigenvalues");
   
  auto alphas = res.tmatrix.diag();
  auto betas = res.tmatrix.offdiag();
  betas.push_back(res.beta);
  auto eigenvalues = res.eigenvalues;
  dumper["Alphas"] << alphas;
  dumper["Betas"] << betas;
  dumper["Eigenvalues"] << eigenvalues;
  dumper.dump();

  lg.out(1, "E0: {}\n", eigenvalues(0));
  
  MPI_Finalize();
  return EXIT_SUCCESS;
}
