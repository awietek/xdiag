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
  bool measure_kinetic = false;

  parse_cmdline(outfile, latticefile, couplingfile, nup, ndown,
		precision, neigenvalue, iters, verbosity, seed,
		measure_kinetic, argc, argv);

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
  double t1 = MPI_Wtime();
  auto H = HubbardModelMPI<double>(bondlist, couplings, qn);
  double t2 = MPI_Wtime();
  lg.out(1, "done. time: {} secs\n", t2-t1); 
  lg.out(1, "dim: {}\n", FormatWithCommas(H.dim())); 
  
  // Define multiplication function
  int iter = 0;
  auto multiply_H = 
    [&H, &iter, verbosity](const VectorMPI<double>& v, VectorMPI<double>& w) 
    {
      bool verbose = (iter==0) && (verbosity > 0);
      double t1 = MPI_Wtime();
      H.apply_hamiltonian(v, w, verbose);
      double t2 = MPI_Wtime();
      lg.out(2, "iter: {}, time: {} secs\n", iter, t2-t1); 
      ++iter;
    };

  // Create normal distributed random start state
  VectorMPI<double> startstate(H.local_dim());
  normal_dist_t<double> dist(0., 1.);
  normal_gen_t<double> gen(dist, seed);
  Random(startstate, gen, true);
  Normalize(startstate);  

  if (measure_kinetic)
    {
      // Reset iters to reasonable size for small model dimension
      iters = std::min(H.dim() / 8 + 5, (unsigned long)iters);
      
      // Define fixed number of steps convergence criterion
      auto converged = 
	[iters](const lila::Tmatrix<double>& tmat, double beta) {
	  (void)beta;
	  return LanczosConvergedFixed(tmat, iters);
	};

      // Define trivial linear combination to get all Lanczos vectors
      std::vector<Vector<double>> linear_combinations;
      for (int i=0; i< iters; ++i)
	{
	  auto lin = Zeros<double>(iters);
	  lin(i) = 1.;
	  linear_combinations.push_back(lin);
	}
      
      // First Lanczos run with random startstate
      auto v0 = startstate;
      auto start_res = Lanczos(multiply_H, v0, converged, 
			       linear_combinations);
      auto alphas_v = start_res.tmatrix.diag();
      auto betas_v = start_res.tmatrix.offdiag();
      betas_v.push_back(start_res.beta);
      auto eigenvalues_v = start_res.eigenvalues;
      auto& vs = start_res.vectors; 
      dumper["Alphas"] << alphas_v;
      dumper["Betas"] << betas_v;
      dumper["Eigenvalues"] << eigenvalues_v;
      dumper["Dimension"] << H.dim();

      Couplings kin_couplings;
      kin_couplings["T"] = 1.;
      kin_couplings["U"] = 0.;
      kin_couplings["C"] = 0.;
      auto T = HubbardModelMPI<double,uint32>(bondlist, kin_couplings, qn);
  
      Matrix<double> vs_T_vs(iters, iters);
      for (int i=0; i<iters; ++i)
	{
	  lg.out(1, "computing kinetic matrix line {}\n", i);
	  auto tmp = vs[i];
	  T.apply_hamiltonian(vs[i], tmp);
	  for (int j=0; j<iters; ++j)
	    vs_T_vs(i,j) = Dot(tmp, vs[j]);
	}
      dumper["VsTVs"] << vs_T_vs;
      dumper.dump();

    }
  else
    {
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
      dumper["Dimension"] << H.dim();
      dumper.dump();

      lg.out(1, "E0: {}\n", eigenvalues(0));
    }
  
  MPI_Finalize();
  return EXIT_SUCCESS;
}
