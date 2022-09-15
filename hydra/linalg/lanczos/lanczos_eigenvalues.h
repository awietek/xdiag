#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/linalg/lanczos/lanczos_convergence.h>
#include <hydra/linalg/lanczos/lanczos_generic.h>
#include <hydra/linalg/lanczos/tmatrix.h>

#include <hydra/blocks/blocks.h>

#ifdef HYDRA_ENABLE_MPI
#include <hydra/parallel/mpi/dot_mpi.h>
#include <hydra/utils/timing_mpi.h>
#else
#include <hydra/utils/timing.h>
#endif

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

// Lanczos which  overwrites starting vector v0
template <class coeff_t, class Block>
Tmatrix LanczosEigenvaluesInplace(
    BondList const &bonds, Couplings const &couplings, Block const &block,
    arma::Col<coeff_t> &v0, int num_eigenvalue = 0, double precision = 1e-12,
    int max_iterations = 1000, double deflation_tol = 1e-7) {

  // MPI Lanczos
#ifdef HYDRA_ENABLE_MPI
  if constexpr (mpi::is_mpi_block<Block>) {
    int iter = 1;
    auto mult = [&iter, &bonds, &couplings, &block](arma::Col<coeff_t> const &v,
                                                    arma::Col<coeff_t> &w) {
      auto ta = rightnow_mpi();
      Apply(bonds, couplings, block, v, block, w);
      Log(1, "Lanczos iteration {}", iter);
      timing_mpi(ta, rightnow_mpi(), "MVM", 1);
      ++iter;
    };

    auto dot_mpi = [](arma::Col<coeff_t> const &v,
                      arma::Col<coeff_t> const &w) -> coeff_t {
      return DotMPI(v, w);
    };

    auto converged = [num_eigenvalue, precision](Tmatrix const &tmat) -> bool {
      return ConvergedEigenvalues(tmat, num_eigenvalue, precision);
    };

    auto t0 = rightnow_mpi();
    auto [tmat, vectors] =
        LanczosGeneric(mult, v0, dot_mpi, converged, arma::Mat<coeff_t>(),
                       max_iterations, deflation_tol);
    (void)vectors;
    timing_mpi(t0, rightnow_mpi(), "Lanczos time", 1);
    return tmat;
  }

  // Serial Lanczos
  else {
#endif
    int iter = 1;
    auto mult = [&iter, &bonds, &couplings, &block](arma::Col<coeff_t> const &v,
                                                    arma::Col<coeff_t> &w) {
      auto ta = rightnow();
      Apply(bonds, couplings, block, v, block, w);
      Log(1, "Lanczos iteration {}", iter);
      timing(ta, rightnow(), "MVM", 1);
      ++iter;
    };

    auto dot = [](arma::Col<coeff_t> const &v,
                  arma::Col<coeff_t> const &w) -> coeff_t {
      return arma::dot(v, w);
    };

    auto converged = [num_eigenvalue, precision](Tmatrix const &tmat) -> bool {
      return ConvergedEigenvalues(tmat, num_eigenvalue, precision);
    };

    auto t0 = rightnow();
    auto [tmat, vectors] =
        LanczosGeneric(mult, v0, dot, converged, arma::Mat<coeff_t>(),
                       max_iterations, deflation_tol);
    (void)vectors;
    timing(t0, rightnow(), "Lanczos time", 1);
    return tmat;
#ifdef HYDRA_ENABLE_MPI
  }
#endif
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
