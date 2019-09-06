#ifndef HYDRA_APPLICATIONS_HUBBARDOPTICAL_ITERATIONSMPI_H_
#define HYDRA_APPLICATIONS_HUBBARDOPTICAL_ITERATIONSMPI_H_

#include <lila/allmpi.h>
#include <hydra/allmpi.h>

#include <mpi.h>

namespace hydra {

struct dyn_lanczos_result_t
{
  lila::Vector<double> alphas;
  lila::Vector<double> betas;
  lila::Vector<double> eigenvalues;
  double dyn_weight;
};


dyn_lanczos_result_t hubbard_dynamical_iterations_lanczos_mpi
(models::HubbardModelMPI<double, uint32>& model,
 const lila::VectorMPI<double>& groundstate, 
 int site1, int site2, std::string fermiontype, int dyniters,
 double precision, int verbosity)
{
  using hydra::models::HubbardModelMPI;
  using namespace lila;

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  
  if ((verbosity >= 1) && (mpi_rank == 0))
    {
      printf("Computing Greens function: type: %s, s1: %d, s2: %d\n",
	     fermiontype.c_str(), site1, site2);
      printf("Applying creation/annihilation operator (%s)...\n",
	     fermiontype.c_str());
      printf("dim before: %ld\n", groundstate.size_global());
    }

  double t1 = MPI_Wtime();
  VectorMPI<double> dyn_start_state1; 
  auto qn_after = 
    model.apply_fermion(groundstate, dyn_start_state1, fermiontype, site1);
  VectorMPI<double> dyn_start_state2; 
  qn_after = 
    model.apply_fermion(groundstate, dyn_start_state2, fermiontype, site2);
  VectorMPI<double> dyn_start_state = dyn_start_state1 + dyn_start_state2;

  dyn_start_state1.clear();
  dyn_start_state2.clear();
  dyn_start_state1.shrink_to_fit();
  dyn_start_state2.shrink_to_fit();

  double dyn_weight = pow(Norm(dyn_start_state), 2);
  double t2 = MPI_Wtime();
  if ((verbosity >= 1) && (mpi_rank == 0))
    {
      printf("dim after: %ld\n", dyn_start_state.size_global());  
      printf("time fermion: %3.4f\n", t2-t1); 
      printf("Creating Hubbard model for n_upspins=%d, n_downspins=%d...\n", 
	     qn_after.n_upspins, qn_after.n_downspins);
    }
  t1 = MPI_Wtime();
  auto model_dyn = model;
  model_dyn.set_qn(qn_after);
  t2 = MPI_Wtime();
  if ((verbosity >= 1) && (mpi_rank == 0))
    { 
      printf("time init dyn: %3.4f\n", t2-t1); 
      printf("Starting dynamical Lanczos procedure ...\n");
    }


  auto multiply_dyn = [&model_dyn, &mpi_rank, &verbosity]
    (const VectorMPI<double>& v, VectorMPI<double>& w) {
    static int iter=0;
    double t1 = MPI_Wtime();
    bool verbose = ((iter == 0) && verbosity >=1);
    model_dyn.apply_hamiltonian(v, w, verbose);
    double t2 = MPI_Wtime();
    if ((verbosity >= 2) && (mpi_rank == 0)) 
      printf("dyniter: %d, time MVM: %3.4f\n", iter, t2-t1); 
    ++iter;
  };

  
  int num_eigenvalue = 1;
  // int random_seed = 0;
  // uint64 dim_dyn = model_dyn.local_dim();
  
  // auto lzs_dyn = Lanczos<double, decltype(multiply_dyn), VectorMPI<double>>
  //   (dim_dyn, random_seed, dyniters, precision, num_eigenvalue, multiply_dyn);
  // lzs_dyn.set_init_state(dyn_start_state);
  // Vector<double> dyn_eigs = lzs_dyn.eigenvalues();

  auto lczs_res = LanczosEigenvalues(multiply_dyn, dyn_start_state, precision, 
				     num_eigenvalue, "Eigenvalues");
  
  
  dyn_lanczos_result_t res;
  // res.alphas = lzs_dyn.alphas();
  // res.betas = lzs_dyn.betas();
  // res.eigenvalues = dyn_eigs;
  // res.dyn_weight = dyn_weight;
  res.alphas = lczs_res.tmatrix.diag();
  res.betas = lczs_res.tmatrix.offdiag();
  res.eigenvalues = lczs_res.eigenvalues;
  res.dyn_weight = dyn_weight;

  return res;
}
  
struct dynamical_iterations_bandlanczos_return_t
{
  lila::Matrix<double> tmatrix;
  lila::Matrix<double> overlaps;
};

dynamical_iterations_bandlanczos_return_t
hubbard_dynamical_iterations_bandlanczos_mpi
(models::HubbardModelMPI<double, uint32>& model,
 const lila::VectorMPI<double>& groundstate, 
 std::vector<int> sites, std::string fermiontype, int dyniters,
 double precision, int verbosity, double deflationtol)
{
  using hydra::models::HubbardModelMPI;
  using namespace lila;

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  
  if ((verbosity >= 1) && (mpi_rank == 0))
    {
      printf("Computing Greens function (BandLanczos): type: %s, \n",
      	fermiontype.c_str());
      printf("Applying creation/annihilation operator (%s)...\n",
      	fermiontype.c_str());
      printf("dim before: %ld\n", groundstate.size_global());
    }

  double t1 = MPI_Wtime();
  std::vector<VectorMPI<double>> dyn_start_states;
  assert(sites.size() > 0);
  int p = sites.size();
  dyn_start_states.resize(p);
  hilbertspaces::hubbard_qn qn_after;
  for (int i=0; i<p; ++i)
    qn_after = model.apply_fermion(groundstate, dyn_start_states[i], fermiontype, sites[i]);
  double t2 = MPI_Wtime();

  if ((verbosity >= 1) && (mpi_rank == 0))
    {
      printf("dim after: %ld\n", dyn_start_states[0].size_global());  
      printf("time fermion: %3.4f\n", t2-t1); 
    }

  // Compute overlap of start states with Lanczos vectors (orthonormal basis thereof)
  auto ortho = lila::gramschmidt(dyn_start_states);
  Matrix<double> overlaps(p, p);
  for (int i=0; i<p; ++i)
    for (int j=0; j<p; ++j)
      {
	// double sdot = Dot(dyn_start_states[i], dyn_start_states[j]);
	// double odot = Dot(ortho[i], ortho[j]);
	// if (mpi_rank == 0) printf("sdot, i: %d, j: %d, %g, \n", i,j,sdot);
	// if (mpi_rank == 0) printf("odot, i: %d, j: %d, %g, \n", i,j,odot);
	// std::cout << std::flush;
	if (i==j) assert(lila::close(Dot(ortho[i], ortho[j]), 1.));
	else assert(lila::close(Dot(ortho[i], ortho[j]), 0.));
	overlaps(i, j) = Dot(ortho[i], dyn_start_states[j]);
      }
  ortho.clear();
  ortho.shrink_to_fit();

  if ((verbosity >= 1) && (mpi_rank == 0)) 
    printf("Creating Hubbard model for n_upspins=%d, n_downspins=%d...\n", 
	   qn_after.n_upspins, qn_after.n_downspins);
  t1 = MPI_Wtime();
  auto model_dyn = model;
  model_dyn.set_qn(qn_after);
  t2 = MPI_Wtime();
  if ((verbosity >= 1) && (mpi_rank == 0))
    {
      printf("time init dyn: %3.4f\n", t2-t1); 
      printf("Starting dynamical Lanczos procedure ...\n");
    }

  auto multiply_dyn = [&model_dyn, &mpi_rank, &verbosity]
    (const VectorMPI<double>& v, VectorMPI<double>& w) {
    static int iter=0;
    double t1 = MPI_Wtime();
    bool verbose = ((iter == 0) && verbosity >=1);
    model_dyn.apply_hamiltonian(v, w, verbose);
    double t2 = MPI_Wtime();
    if ((verbosity >= 2) && (mpi_rank == 0)) 
      printf("dyniter: %d, time MVM: %3.4f\n", iter, t2-t1); 
    ++iter;
  };

  uint64 dim_dyn = model_dyn.local_dim();
  int random_seed = 0;
  int num_eigenvalue = 1;
  auto lzs_dyn = BandLanczos<double, decltype(multiply_dyn), VectorMPI<double>>
    (dim_dyn, random_seed, dyniters, precision, num_eigenvalue, multiply_dyn, p);
  lzs_dyn.set_init_states(dyn_start_states);

  auto res = lzs_dyn.eigenvalues();
  dynamical_iterations_bandlanczos_return_t ret;
  ret.tmatrix = lzs_dyn.tmatrix();
  ret.overlaps = overlaps;
  return ret;
}


}
#endif
