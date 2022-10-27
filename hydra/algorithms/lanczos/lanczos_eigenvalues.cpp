#include "lanczos_eigenvalues.h"

#include <hydra/algebra/algebra.h>
#include <hydra/algorithms/lanczos/lanczos.h>
#include <hydra/algorithms/lanczos/lanczos_convergence.h>
#include <hydra/states/random_state.h>
#include <hydra/utils/timing.h>

namespace hydra {

template <class coeff_t>
Tmatrix lanczos_eigenvalues_inplace(BondList const &bonds,
                                    State<coeff_t> &state_0, int num_eigenvalue,
                                    double precision, int max_iterations,
                                    double deflation_tol) {

  auto const &block = state_0.block();
  auto &v0 = state_0.vector();

  int iter = 1;
  auto mult = [&iter, &bonds, &block](arma::Col<coeff_t> const &v,
                                      arma::Col<coeff_t> &w) {
    auto ta = rightnow();
    apply(bonds, block, v, block, w);
    Log(1, "Lanczos iteration {}", iter);
    timing(ta, rightnow(), "MVM", 1);
    ++iter;
  };

  auto converged = [num_eigenvalue, precision](Tmatrix const &tmat) -> bool {
    return converged_eigenvalues(tmat, num_eigenvalue, precision);
  };

  auto t0 = rightnow();
  auto tmat = lanczos(mult, v0, converged, max_iterations, deflation_tol);
  timing(t0, rightnow(), "Lanczos time", 1);
  return tmat;
}

template Tmatrix lanczos_eigenvalues_inplace(BondList const &, State<double> &,
                                             int, double, int, double);
template Tmatrix lanczos_eigenvalues_inplace(BondList const &, State<complex> &,
                                             int, double, int, double);

template <class coeff_t>
Tmatrix lanczos_eigenvalues(BondList const &bonds, State<coeff_t> state_0,
                            int num_eigenvalue, double precision,
                            int max_iterations, double deflation_tol) {
  return lanczos_eigenvalues_inplace(bonds, state_0, num_eigenvalue, precision,
                                     max_iterations, deflation_tol);
}

template Tmatrix lanczos_eigenvalues(BondList const &, State<double>, int,
                                     double, int, double);
template Tmatrix lanczos_eigenvalues(BondList const &, State<complex>, int,
                                     double, int, double);

// Lanczos which does not overwrite v0
template <class coeff_t>
Tmatrix lanczos_eigenvalues(BondList const &bonds, Block const &block,
                            int num_eigenvalue, double precision,
                            int max_iterations, double deflation_tol) {
  auto rstate = RandomState();
  auto state_0 = State<coeff_t>(block, rstate);
  return lanczos_eigenvalues_inplace(bonds, state_0, num_eigenvalue, precision,
                                     max_iterations, deflation_tol);
}

template Tmatrix lanczos_eigenvalues<double>(BondList const &, Block const &,
                                             int, double, int, double);
template Tmatrix lanczos_eigenvalues<complex>(BondList const &, Block const &,
                                              int, double, int, double);

Tmatrix lanczos_eigenvalues_real(BondList const &bonds, Block const &block,
                                 int num_eigenvalue, double precision,
                                 uint64_t seed, int max_iterations,
                                 double deflation_tol) {
  auto rstate = RandomState(seed);
  auto state_0 = StateReal(block, rstate);
  return lanczos_eigenvalues_inplace(bonds, state_0, num_eigenvalue, precision,
                                     max_iterations, deflation_tol);
}

Tmatrix lanczos_eigenvalues_cplx(BondList const &bonds, Block const &block,
                                 int num_eigenvalue, double precision,
                                 uint64_t seed, int max_iterations,
                                 double deflation_tol) {
  auto rstate = RandomState(seed);
  auto state_0 = StateCplx(block, rstate);
  return lanczos_eigenvalues_inplace(bonds, state_0, num_eigenvalue, precision,
                                     max_iterations, deflation_tol);
}

} // namespace hydra
