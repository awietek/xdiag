#include <cstdlib>

#include <mpi.h>

#include <lila/allmpi.h>
#include <hydra/allmpi.h>
#include <lime/all.h>

lila::LoggerMPI lg;

#include "tjed.options.h"


template <class coeff_t>
void run_real_complex(std::string real_complex,
		      hydra::all::BondList bondlist,
		      hydra::all::Couplings couplings,
		      hydra::all::hubbard_qn qn,
		      double precision, int neigenvalue,
		      int iters, int verbosity, int seed,
		      std::string outfile)
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  
  auto dumper = lime::MeasurementsH5((mpi_rank == 0) ? outfile : "");

  lg.out(1, "Creating {} t-J model for n_upspins={}, n_downspins={}...\n",
	 real_complex, qn.n_upspins, qn.n_downspins);
  lg.out(1, "Using {} MPI tasks\n", mpi_size);

  double t1 = MPI_Wtime();
  auto H = TJModelMPI<coeff_t>(bondlist, couplings, qn);
  double t2 = MPI_Wtime();
  lg.out(1, "done. time: {} secs\n", t2-t1); 
  lg.out(1, "dim: {}\n", FormatWithCommas(H.dim())); 


  // Define multiplication function
  int iter = 0;
  auto multiply_H = 
    [&H, &iter, verbosity](const VectorMPI<coeff_t>& v, VectorMPI<coeff_t>& w) 
    {
      bool verbose = (iter==0) && (verbosity > 0);
      double t1 = MPI_Wtime();
      H.apply_hamiltonian(v, w, verbose);
      double t2 = MPI_Wtime();
      lg.out(2, "iter: {}, time: {} secs\n", iter, t2-t1); 
      ++iter;
    };

  // Create normal distributed random start state
  VectorMPI<coeff_t> startstate(H.local_dim());
  normal_dist_t<coeff_t> dist(0., 1.);
  normal_gen_t<coeff_t> gen(dist, seed);
  Random(startstate, gen, true);
  Normalize(startstate);  

  // Run Lanczos
  auto res = LanczosEigenvalues(multiply_H, startstate, precision,
				neigenvalue, "Eigenvalues", iters);
  if (res.exhausted) 
    lg.out("Warning: Lanczos sequence exhausted after {} steps\n", res.eigenvalues.size());
  if (!res.converged) 
    lg.out("Warning: Lanczos sequence not converged in {} steps\n", iters);
   
  auto alphas = res.tmatrix.diag();
  auto betas = res.tmatrix.offdiag();
  betas.push_back(res.beta);
  auto eigenvalues = res.eigenvalues;
  dumper["Alphas"] << alphas;
  dumper["Betas"] << betas;
  dumper["Eigenvalues"] << eigenvalues;
  dumper["Dimension"] << H.dim();
  dumper.dump();

  lg.out(1, "E0: {}\n", eigenvalues(0));
  
}

int main(int argc, char* argv[])
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;

  MPI_Init(&argc, &argv); 

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


  if (couplings.all_real())
    run_real_complex<double>(std::string("REAL"), bondlist,
			     couplings, qn, precision, neigenvalue, iters,
			     verbosity, seed, outfile);
  else
    run_real_complex<lila::complex>(std::string("COMPLEX"), bondlist,
				    couplings, qn, precision, neigenvalue,
				    iters, verbosity, seed, outfile);
  
  MPI_Finalize();
  return EXIT_SUCCESS;
}
