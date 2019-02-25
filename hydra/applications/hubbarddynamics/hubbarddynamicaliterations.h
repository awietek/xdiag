#ifndef HYDRA_APPLICATIONS_HUBBARDDYNAMICS_ITERATIONS_H_
#define HYDRA_APPLICATIONS_HUBBARDDYNAMICS_ITERATIONS_H_
#include <chrono>

#include <lila/all.h>
#include <hydra/all.h>

namespace hydra {

struct dyn_lanczos_result_t
{
  lila::Vector<double> alphas;
  lila::Vector<double> betas;
  lila::Vector<double> eigenvalues;
  double dyn_weight;
};


dyn_lanczos_result_t hubbard_dynamical_iterations_lanczos
(models::HubbardModel& model,
 const lila::Vector<double>& groundstate, 
 int site1, int site2, std::string fermiontype, int dyniters)
{
  using hydra::models::HubbardModel;
  using namespace lila;
  using Clock = std::chrono::high_resolution_clock;
  using secs = std::chrono::duration<float>;

  printf("Computing Greens function: type: %s, s1: %d, s2: %d\n",
			    fermiontype.c_str(), site1, site2);

  printf("Applying creation/annihilation operator (%s)...\n",
			    fermiontype.c_str());
  auto t1 = Clock::now();
  printf("dim before: %d\n", (int)groundstate.size());
  Vector<double> dyn_start_state1; 
  auto qn_after = 
    model.apply_fermion(groundstate, dyn_start_state1, fermiontype, site1);
  Vector<double> dyn_start_state2; 
  qn_after = 
    model.apply_fermion(groundstate, dyn_start_state2, fermiontype, site2);
  Vector<double> dyn_start_state = dyn_start_state1 + dyn_start_state2;

  dyn_start_state1.clear();
  dyn_start_state2.clear();
  dyn_start_state1.shrink_to_fit();
  dyn_start_state2.shrink_to_fit();

  double dyn_weight = pow(Norm(dyn_start_state), 2);
  printf("dim after: %d\n", (int)dyn_start_state.size());  
  auto t2 = Clock::now();
  printf("time fermion: %3.4f\n", secs(t2-t1).count()); 

  
  printf("Creating Hubbard model for n_upspins=%d, n_downspins=%d...\n", 
	 qn_after.n_upspins, qn_after.n_downspins);
  t1 = Clock::now();
  auto model_dyn = model;
  model_dyn.set_qn(qn_after);
  t2 = Clock::now();
  printf("time init dyn: %3.4f\n", secs(t2-t1).count()); 


  printf("Starting dynamical Lanczos procedure ...\n");
  double precision = -1;
  int max_iterations = dyniters;

  auto multiply_dyn = [&model_dyn]
    (const Vector<double>& v, Vector<double>& w) {
    static int iter=0;
    auto t1 = Clock::now();
    model_dyn.apply_hamiltonian(v, w);
    auto t2 = Clock::now();
    printf("dyniter: %d, time MVM: %3.4f\n", iter, secs(t2-t1).count()); 
    ++iter;
  };

  uint64 dim_dyn = model_dyn.dim();
  int random_seed = 0;
  int num_eigenvalue = 1;
  auto lzs_dyn = Lanczos<double, decltype(multiply_dyn), Vector<double>>
    (dim_dyn, random_seed, max_iterations, precision, num_eigenvalue, multiply_dyn);
  lzs_dyn.set_init_state(dyn_start_state);

  Vector<double> dyn_eigs = lzs_dyn.eigenvalues();
  
  dyn_lanczos_result_t res;
  res.alphas = lzs_dyn.alphas();
  res.betas = lzs_dyn.betas();
  res.eigenvalues = dyn_eigs;
  res.dyn_weight = dyn_weight;
  return res;
}
 

struct dynamical_iterations_bandlanczos_return_t
{
  lila::Matrix<double> tmatrix;
  lila::Matrix<double> overlaps;
};

dynamical_iterations_bandlanczos_return_t
hubbard_dynamical_iterations_bandlanczos
(models::HubbardModel& model, const lila::Vector<double>& groundstate, 
 std::vector<int> sites, std::string fermiontype, int dyniters)
{
  using hydra::models::HubbardModel;
  using namespace lila;
  using Clock = std::chrono::high_resolution_clock;
  using secs = std::chrono::duration<float>;

  printf("Computing Greens function (BandLanczos): type: %s, \n",
			    fermiontype.c_str());

  printf("Applying creation/annihilation operator (%s)...\n",
			    fermiontype.c_str());
  auto t1 = Clock::now();
  printf("dim before: %d\n", (int)groundstate.size());
  std::vector<Vector<double>> dyn_start_states;
  assert(sites.size() > 0);
  int p = sites.size();
  dyn_start_states.resize(p);
  hilbertspaces::hubbard_qn qn_after;
  for (int i=0; i<p; ++i)
    qn_after = model.apply_fermion(groundstate, dyn_start_states[i], fermiontype, sites[i]);

  printf("dim after: %d\n", (int)dyn_start_states[0].size());  
  auto t2 = Clock::now();
  printf("time fermion: %3.4f\n", secs(t2-t1).count()); 

  // Compute overlap of start states with Lanczos vectors (orthonormal basis thereof)
  auto ortho = lila::gramschmidt(dyn_start_states);
  Matrix<double> overlaps(p, p);
  for (int i=0; i<p; ++i)
    for (int j=0; j<p; ++j)
      overlaps(i, j) = Dot(ortho[i], dyn_start_states[j]);
  ortho.clear();
  ortho.shrink_to_fit();

  printf("Creating Hubbard model for n_upspins=%d, n_downspins=%d...\n", 
	 qn_after.n_upspins, qn_after.n_downspins);
  t1 = Clock::now();
  auto model_dyn = model;
  model_dyn.set_qn(qn_after);
  t2 = Clock::now();
  printf("time init dyn: %3.4f s\n", secs(t2-t1).count()); 


  printf("Starting dynamical Lanczos procedure ...\n");
  double precision = -1;
  int max_iterations = dyniters;

  auto multiply_dyn = [&model_dyn]
    (const Vector<double>& v, Vector<double>& w) {
    static int iter=0;
    auto t1 = Clock::now();
    model_dyn.apply_hamiltonian(v, w);
    auto t2 = Clock::now();
    printf("dyniter: %d, time MVM: %3.4f s\n", iter, secs(t2-t1).count()); 
    ++iter;
  };

  uint64 dim_dyn = model_dyn.dim();
  int random_seed = 0;
  int num_eigenvalue = 1;
  auto lzs_dyn = BandLanczos<double, decltype(multiply_dyn), Vector<double>>
    (dim_dyn, random_seed, max_iterations, precision, num_eigenvalue, multiply_dyn, p);
  lzs_dyn.set_init_states(dyn_start_states);

  auto res = lzs_dyn.eigenvalues();
  dynamical_iterations_bandlanczos_return_t ret;
  ret.tmatrix = lzs_dyn.tmatrix();
  ret.overlaps = overlaps;
  return ret;
}


}
#endif
