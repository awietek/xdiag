#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/algorithms/lanczos/lanczos.h>
#include <hydra/algorithms/lanczos/lanczos_convergence.h>
#include <hydra/algorithms/lanczos/tmatrix.h>
#include <hydra/blocks/blocks.h>
#include <hydra/operators/bondlist.h>
#include <hydra/utils/timing.h>

namespace hydra {

// Lanczos which  overwrites starting vector v0
template <class coeff_t, class Block>
Tmatrix
lanczos_eigenvalues_inplace(BondList const &bonds, Block const &block,
                            arma::Col<coeff_t> &v0, int num_eigenvalue = 0,
                            double precision = 1e-12, int max_iterations = 1000,
                            double deflation_tol = 1e-7) {

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

// Lanczos which does not overwrite v0
template <class coeff_t, class Block>
Tmatrix lanczos_eigenvalues(BondList const &bonds, Block const &block,
                            arma::Col<coeff_t> v0, int num_eigenvalue = 0,
                            double precision = 1e-12, int max_iterations = 1000,
                            double deflation_tol = 1e-7) {
  return lanczos_eigenvalues_inplace(bonds, block, v0, num_eigenvalue,
                                     precision, max_iterations, deflation_tol);
}

template <class Block>
Tmatrix lanczos_eigenvalues_real(BondList const &bonds, Block const &block,
                                 int num_eigenvalue = 0,
                                 double precision = 1e-12, int seed = 42,
                                 int max_iterations = 1000,
                                 double deflation_tol = 1e-7) {
  auto v0 = random_state_real(block, seed);
  return lanczos_eigenvalues_inplace(bonds, block, v0.vector(), num_eigenvalue,
                                     precision, max_iterations, deflation_tol);
}

template <class Block>
Tmatrix lanczos_eigenvalues_cplx(BondList const &bonds, Block const &block,
                                 int num_eigenvalue = 0,
                                 double precision = 1e-12, int seed = 42,
                                 int max_iterations = 1000,
                                 double deflation_tol = 1e-7) {
  auto v0 = random_state_cplx(block, seed);
  return lanczos_eigenvalues_inplace(bonds, block, v0.vector(), num_eigenvalue,
                                     precision, max_iterations, deflation_tol);
}

} // namespace hydra
