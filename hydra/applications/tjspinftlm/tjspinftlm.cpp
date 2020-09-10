#include <cstdlib>
#include <vector>

#include <mpi.h>
#include <hydra/allmpi.h>
#include <lila/allmpi.h>
#include <lime/all.h>

#include "tjspinftlm.options.h"

lila::LoggerMPI lg;

template <class coeff_t>
void run_real_complex(std::string real_complex,
    hydra::all::BondList bondlist,
    hydra::all::Couplings couplings,
    hydra::all::BondList corr_bondlist,
    hydra::all::hubbard_qn qn,
    int iters, 
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

  lg.out(1, "Creating tJ model for n_upspins={}, n_downspins={}...\n",
	     qn.n_upspins, qn.n_downspins);
  double t1 = MPI_Wtime();
  auto H = TJModelMPI<coeff_t>(bondlist, couplings, qn);
  double t2 = MPI_Wtime();
  lg.out(1, "time init: {} secs\n", t2-t1); 


  // Define multiplication function
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
  
  // Create normal distributed random start state
  int random_seed = seed + 1234567*mpi_rank;
  VectorMPI<coeff_t> startstate(H.local_dim());
  normal_dist_t<coeff_t> dist(0., 1.);
  normal_gen_t<coeff_t> gen(dist, random_seed);
  Random(startstate, gen, true);
  Normalize(startstate);



  // Reset iters to reasonable size for small model dimension
  iters = std::min(H.dim() / 8 + 5, (unsigned long)iters);

  // Define fixed number of steps convergence criterion
  auto converged = 
    [iters](const lila::Tmatrix<double>& tmat, double beta) {
    (void)beta;
    return LanczosConvergedFixed(tmat, iters);
  };


  // Define trivial linear combination to get all Lanczos vectors
  std::vector<Vector<coeff_t>> linear_combinations;
  for (int i=0; i< iters; ++i)
    {
      auto lin = Zeros<coeff_t>(iters);
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
  if (mpi_rank == 0) {
  file["AlphasV"] = alphas_v;
  file["BetasV"] = betas_v;
  file["EigenvaluesV"] = eigenvalues_v;
  }
  
  // Compute overlap of Lanczos vectors with start vectors
   
  Vector<coeff_t> r_dot_vs(iters);
  for (int i=0; i<iters; ++i)
    {
      r_dot_vs(i) = Dot(startstate, vs[i]);
    }
  if (mpi_rank == 0) {
  file["RDotVs"] = r_dot_vs;
  }
  
  // Calculate eigenvalues/vectors of T matrix
  
  auto tmatrix_results = Eigen(start_res.tmatrix);
  auto tmatrix_eigenvectors = tmatrix_results.eigenvectors;

  // Assemble ground state
  
  VectorMPI<coeff_t> groundstate(H.local_dim());
  Zeros(groundstate.vector_local());
  iters = alphas_v.size();
  
  for (int i=0; i < iters; i++) {
    auto tmpVec = vs[i];
    coeff_t coeff = tmatrix_eigenvectors(0, i);
    Scale(coeff, tmpVec);
    groundstate += tmpVec;
  }
  groundstate = vs[0];
  
  // Start loop through all bonds
  auto corrs_zs = corr_bondlist.bonds_of_type("SZ");
  for (auto corr : corrs_zs)
  {
    assert(corr.sites().size() == 2);
    int s1 = corr.sites()[0];
    int s2 = corr.sites()[1];
    std::stringstream ss;
    ss << "sz_" << s1 << "_" << s2;
    std::string label = ss.str();

    // Second Lanczos run with current startstate
    iter = 0;
    H.apply_sz(groundstate, v0, s1);
    auto sz_start_res = Lanczos(multiply_H, v0, converged, 
             linear_combinations);
    auto alphas_v_tilde = sz_start_res.tmatrix.diag();
    auto betas_v_tilde = sz_start_res.tmatrix.offdiag();
    betas_v_tilde.push_back(sz_start_res.beta);
    auto eigenvalues_v_tilde = sz_start_res.eigenvalues;
    auto& vs_tilde = sz_start_res.vectors;
    if (mpi_rank == 0)
    {
    file[label + std::string("_AlphasVTilde")] = alphas_v_tilde;
    file[label + std::string("_BetasVTilde")] = betas_v_tilde;
    file[label + std::string("_EigenvaluesVTilde")] = eigenvalues_v_tilde;
    }

    int iters_tilde = alphas_v_tilde.size();


    // Compute overlap of sz Lanczos vectors with sz start vectors
    Vector<coeff_t> vs_tilde_dot_A_r(iters_tilde);
    for (int i=0; i<iters_tilde; ++i)
      {
        vs_tilde_dot_A_r(i) = Dot(vs_tilde[i], groundstate);
      }
    if (mpi_rank == 0) {
    file[label + std::string("_VsTildeDotAR")] = vs_tilde_dot_A_r;
    }


    // Compute sz operator matrix
    Matrix<coeff_t> vs_A_vs_tilde(iters, iters_tilde);
    for (int i=0; i<iters; ++i)
      {
        auto tmp = vs[i];
        H.apply_sz(vs[i], tmp, s2);
        for (int j=0; j<iters_tilde; ++j)
    vs_A_vs_tilde(i,j) = Dot(tmp, vs_tilde[j]);
      }

    
    if (mpi_rank == 0){
    file[label + std::string("_VsAVsTilde")] = vs_A_vs_tilde;
    }
  }

}
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
  int nup = -1;
  int ndown = -1;
  int iters = 1000;
  int verbosity = 1;
  int seed = 1;
  parse_cmdline(outfile, latticefile, couplingfile, corrfile, 
		nup, ndown, iters, verbosity, seed, argc, argv);
  lg.set_verbosity(verbosity);  
  lg.out(1, "Using {} MPI tasks\n", mpi_size);

  check_if_files_exists({latticefile, couplingfile, corrfile});

  // Parse bondlist, couplings, and corr_bondlist from file
  BondList bondlist = read_bondlist(latticefile);
  Couplings couplings = read_couplings(couplingfile);
  BondList corr_bondlist = read_bondlist(corrfile);

  // Create infrastructure for tJ model
  int n_sites = bondlist.n_sites();
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};

  if (couplings.all_real())
    run_real_complex<double>(std::string("REAL"), bondlist,
        couplings, corr_bondlist, qn, iters,
        verbosity, seed, outfile);
  else
    run_real_complex<lila::complex>(std::string("COMPLEX"), bondlist,
        couplings, corr_bondlist, qn, iters,
         verbosity, seed, outfile);

  MPI_Finalize();
  return EXIT_SUCCESS;
}

