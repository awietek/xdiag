#ifndef HYDRA_APPLICATIONS_TJDYNAMICS_ITERATIONSMPI_H_
#define HYDRA_APPLICATIONS_TJDYNAMICS_ITERATIONSMPI_H_

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

template <class coeff_t>
dyn_lanczos_result_t tj_dynamical_iterations_lanczos_mpi
(models::TJModelMPI<coeff_t>& model,
 const lila::VectorMPI<coeff_t>& groundstate, 
 int site1, int site2, std::string operatortype, int dyniters,
 double precision)
{
  using hydra::models::TJModelMPI;
  using namespace lila;

  lg.out(1, "Computing Greens function: type: {}, s1: {}, s2: {}\n",
	 operatortype.c_str(), site1, site2);

  double t1 = MPI_Wtime();
  VectorMPI<coeff_t> dyn_start_state1; 
    model.apply_sz(groundstate, dyn_start_state1, site1);
  VectorMPI<coeff_t> dyn_start_state2; 
    model.apply_sz(groundstate, dyn_start_state2, site2);
  VectorMPI<coeff_t> dyn_start_state = dyn_start_state1 + dyn_start_state2;
/*
  dyn_start_state1.clear();
  dyn_start_state2.clear();
  dyn_start_state1.shrink_to_fit();
  dyn_start_state2.shrink_to_fit();
*/
  double dyn_weight = pow(Norm(dyn_start_state), 2);
  double t2 = MPI_Wtime();
  lg.out(1, "time sz: {}\n", t2-t1); 


  auto model_dyn = model;

  // Define multiplication function
  int dyn_iter = 0;
  auto multiply_dyn = [&model_dyn, &dyn_iter](const VectorMPI<coeff_t>& v, 
					      VectorMPI<coeff_t>& w) 
    {
      double t1 = MPI_Wtime();
      model_dyn.apply_hamiltonian(v, w);
      double t2 = MPI_Wtime();
      lg.out(2, "dyn iter: {}, time MVM: {}\n", dyn_iter, t2-t1); 
      ++dyn_iter;
    };

  // Run dynamical Lanczos
  lg.out(1, "Starting dynamical Lanczos procedure ...\n");
  int num_eigenvalue = 1;
  auto lczs_res = LanczosEigenvalues(multiply_dyn, dyn_start_state, precision, 
				     num_eigenvalue, "Eigenvalues");
 
  
  dyn_lanczos_result_t res;
  res.alphas = lczs_res.tmatrix.diag();
  res.betas = lczs_res.tmatrix.offdiag();
  res.eigenvalues = lczs_res.eigenvalues;
  res.dyn_weight = dyn_weight;

  return res;
}
/*  
struct dynamical_iterations_bandlanczos_return_t
{
  lila::Matrix<double> tmatrix;
  lila::Matrix<double> overlaps;
};

template <class coeff_t>
dynamical_iterations_bandlanczos_return_t
tj_dynamical_iterations_bandlanczos_mpi
(models::tJModelMPI<coeff_t>& model,
 const lila::VectorMPI<coeff_t>& groundstate, 
 std::vector<int> sites, std::string operatortype, int dyniters,
 double precision, int verbosity, double deflationtol)
{
  using hydra::models::HubbardModelMPI;
  using namespace lila;

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  
  if ((verbosity >= 1) && (mpi_rank == 0))
    {
      printf("Applying sz operator...\n");
      printf("dim before: %ld\n", groundstate.size_global());
    }

  double t1 = MPI_Wtime();
  std::vector<VectorMPI<coeff_t>> dyn_start_states;
  assert(sites.size() > 0);
  int p = sites.size();
  dyn_start_states.resize(p);
  hilbertspaces::hubbard_qn qn_after;
  for (int i=0; i<p; ++i)
    qn_after = model.apply_sz(groundstate, dyn_start_states[i], sites[i]);
  double t2 = MPI_Wtime();

  if ((verbosity >= 1) && (mpi_rank == 0))
    {
      printf("time sz: %3.4f\n", t2-t1); 
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
*/

//template class tj_dynamical_iterations_lanczos_mpi<double>;
//template class tj_dynamical_iterations_lanczos_mpi<complex>;

}
#endif
