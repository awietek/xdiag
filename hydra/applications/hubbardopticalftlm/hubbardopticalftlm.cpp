#include <cstdlib>
#include <vector>

#include <mpi.h>
#include <hydra/allmpi.h>
#include <lila/allmpi.h>
#include <lime/all.h>

#include "hubbardopticalftlm.options.h"

lila::LoggerMPI lg;

int main(int argc, char* argv[])
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;

  MPI_Init(&argc, &argv); 
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Parse input
  std::string outfile;
  std::string latticefile;
  std::string couplingfile;
  std::string corrfile;
  std::string corrcouplingfile;
  int nup = -1;
  int ndown = -1;
  int iters = 100;
  int verbosity = 1;
  int seed = 1;
  bool kinetic = false;
  parse_cmdline(outfile, latticefile, couplingfile, corrfile, corrcouplingfile,
		nup, ndown, iters, verbosity, seed, kinetic, argc, argv);
  lg.set_verbosity(verbosity);  
  lg.out(1, "Using {} MPI tasks\n", mpi_size);

  auto dumper = lime::MeasurementsH5((mpi_rank == 0) ? outfile : "");
  check_if_files_exists({latticefile, couplingfile, corrfile, corrcouplingfile});

  // Create Hamiltonian
  BondList bondlist = read_bondlist(latticefile);
  Couplings couplings = read_couplings(couplingfile);

  if (mpi_rank==0) HydraPrint(bondlist);

  // Create infrastructure for Hubbard model
  int n_sites = bondlist.n_sites();
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};

  lg.out(1, "Creating Hubbard model for n_upspins={}, n_downspins={}...\n",
	     qn.n_upspins, qn.n_downspins);
  double t1 = MPI_Wtime();
  auto H = HubbardModelMPI<double,uint32>(bondlist, couplings, qn);
  double t2 = MPI_Wtime();
  lg.out(1, "time init: {} secs\n", t2-t1); 


  // Define multiplication function
  int iter = 0;
  auto multiply_H = 
    [&H, &iter](const VectorMPI<double>& v, VectorMPI<double>& w) 
    {
      double t1 = MPI_Wtime();
      H.apply_hamiltonian(v, w, false);
      double t2 = MPI_Wtime();
      lg.out(2, "iter: {}, time MVM: %3.4f\n", iter, t2-t1); 
      ++iter;
    };


  // Create current operator
  BondList curr_bondlist = read_bondlist(corrfile);
  Couplings curr_couplings = read_couplings(corrcouplingfile);  
  
  lg.out(1, "Creating Current for n_upspins={}, n_downspins={}...\n",
	     qn.n_upspins, qn.n_downspins);
  t1 = MPI_Wtime();
  auto C = HubbardModelMPI<double,uint32>(curr_bondlist, curr_couplings, qn);
  t2 = MPI_Wtime();
  lg.out(1, "time curr init: {}\n", t2-t1); 

  
  // Create normal distributed random start state
  int random_seed = seed + 1234567*mpi_rank;
  VectorMPI<double> startstate(H.local_dim());
  normal_dist_t<double> dist(0., 1.);
  normal_gen_t<double> gen(dist, random_seed);
  Random(startstate, gen, true);
  Normalize(startstate);


  // Create random current start state
  VectorMPI<double> currstartstate = startstate;
  C.apply_hamiltonian(startstate, currstartstate);


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
  dumper["AlphasV"] << alphas_v;
  dumper["BetasV"] << betas_v;
  dumper["EigenvaluesV"] << eigenvalues_v;
  dumper.dump();

  // Second Lanczos run with current startstate
  iter = 0;
  v0 = currstartstate;
  auto current_start_res = Lanczos(multiply_H, v0, converged, 
				   linear_combinations);
  auto alphas_v_tilde = current_start_res.tmatrix.diag();
  auto betas_v_tilde = current_start_res.tmatrix.offdiag();
  betas_v_tilde.push_back(current_start_res.beta);
  auto eigenvalues_v_tilde = current_start_res.eigenvalues;
  auto& vs_tilde = current_start_res.vectors;
  dumper["AlphasVTilde"] << alphas_v_tilde;
  dumper["BetasVTilde"] << betas_v_tilde;
  dumper["EigenvaluesVTilde"] << eigenvalues_v_tilde;
  dumper.dump();

  iters = alphas_v.size();
  int iters_tilde = alphas_v_tilde.size();

  // Compute overlap of Lanczos vectors with start vectors
  Vector<double> r_dot_vs(iters);
  for (int i=0; i<iters; ++i)
    {
      r_dot_vs(i) = Dot(startstate, vs[i]);
      lg.out(1, "r_dot_vs({})= {}\n", i, r_dot_vs(i));
    }
  dumper["RDotVs"] << r_dot_vs;
  dumper.dump();

  // Compute overlap of current Lanczos vectors with current start vectors
  Vector<double> vs_tilde_dot_A_r(iters_tilde);
  for (int i=0; i<iters_tilde; ++i)
    {
      vs_tilde_dot_A_r(i) = Dot(vs_tilde[i], currstartstate);
      lg.out(1, "vs_tilde_dot_A_r({})= {}\n", i, vs_tilde_dot_A_r(i));
    }
  dumper["VsTildeDotAR"] << vs_tilde_dot_A_r;
  dumper.dump();  


  // Compute current operator matrix
  Matrix<double> vs_A_vs_tilde(iters, iters_tilde);
  for (int i=0; i<iters; ++i)
    {
      auto tmp = vs[i];
      lg.out(1, "computing matrix line {}\n", i);
      C.apply_hamiltonian(vs[i], tmp);
      for (int j=0; j<iters_tilde; ++j)
	vs_A_vs_tilde(i,j) = Dot(tmp, vs_tilde[j]);
    }
  dumper["VsAVsTilde"] << vs_A_vs_tilde;
  dumper.dump();


  // Compute kinetic energy matrix
  if (kinetic)
    {
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

  MPI_Finalize();
  return EXIT_SUCCESS;
}
