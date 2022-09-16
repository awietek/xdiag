#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/linalg/lanczos/lanczos_convergence.h>
#include <hydra/linalg/lanczos/lanczos.h>
#include <hydra/linalg/lanczos/tmatrix.h>

#include <hydra/blocks/blocks.h>
#include <hydra/utils/timing.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

// Lanczos which  overwrites starting vector v0
template <class coeff_t, class Block>
Tmatrix LanczosEigenvaluesInplace(
    BondList const &bonds, Couplings const &couplings, Block const &block,
    arma::Col<coeff_t> &v0, int num_eigenvalue = 0, double precision = 1e-12,
    int max_iterations = 1000, double deflation_tol = 1e-7) {

  int iter = 1;
  auto mult = [&iter, &bonds, &couplings, &block](arma::Col<coeff_t> const &v,
                                                  arma::Col<coeff_t> &w) {
    auto ta = rightnow();
    Apply(bonds, couplings, block, v, block, w);
    Log(1, "Lanczos iteration {}", iter);
    timing(ta, rightnow(), "MVM", 1);
    ++iter;
  };

  auto converged = [num_eigenvalue, precision](Tmatrix const &tmat) -> bool {
    return ConvergedEigenvalues(tmat, num_eigenvalue, precision);
  };

  auto t0 = rightnow();
  auto tmat = Lanczos(mult, v0, converged, max_iterations, deflation_tol);
  timing(t0, rightnow(), "Lanczos time", 1);
  return tmat;
}

// Lanczos which does not overwrite v0
template <class coeff_t, class Block>
Tmatrix LanczosEigenvalues(BondList const &bonds, Couplings const &couplings,
                           Block const &block, arma::Col<coeff_t> v0,
                           int num_eigenvalue = 0, double precision = 1e-12,
                           int max_iterations = 1000,
                           double deflation_tol = 1e-7) {

  return LanczosEigenvaluesInplace(bonds, couplings, block, v0, num_eigenvalue,
                                   precision, max_iterations, deflation_tol);
}

template <class Block>
Tmatrix LanczosEigenvaluesReal(BondList const &bonds,
                               Couplings const &couplings, Block const &block,
                               int num_eigenvalue = 0, double precision = 1e-12,
                               int seed = 42, int max_iterations = 1000,
                               double deflation_tol = 1e-7) {

  auto v0 = RandomStateReal(block, seed);
  return LanczosEigenvaluesInplace(bonds, couplings, block, v0.vector(),
                                   num_eigenvalue, precision, max_iterations,
                                   deflation_tol);
}

template <class Block>
Tmatrix LanczosEigenvaluesCplx(BondList const &bonds,
                               Couplings const &couplings, Block const &block,
                               int num_eigenvalue = 0, double precision = 1e-12,
                               int seed = 42, int max_iterations = 1000,
                               double deflation_tol = 1e-7) {

  auto v0 = RandomStateCplx(block, seed);
  return LanczosEigenvaluesInplace(bonds, couplings, block, v0.vector(),
                                   num_eigenvalue, precision, max_iterations,
                                   deflation_tol);
}

} // namespace hydra
