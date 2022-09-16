#pragma once

#include <algorithm>

#include "extern/armadillo/armadillo"

#include <hydra/linalg/lanczos/lanczos_convergence.h>
#include <hydra/linalg/lanczos/lanczos.h>
#include <hydra/linalg/lanczos/lanczos_build.h>
#include <hydra/linalg/lanczos/tmatrix.h>

#include <hydra/blocks/blocks.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <hydra/random/hash_functions.h>
#include <hydra/random/hashes.h>
#include <hydra/random/random_utils.h>

namespace hydra {

// Implementation with a function setting the starting vector
template <class coeff_t, class Block, class set_v0_f>
std::pair<Tmatrix, arma::Col<coeff_t>>
LanczosEigenvector(BondList const &bonds, Couplings const &couplings,
                   Block const &block, set_v0_f set_v0, int num_eigenvector = 0,
                   double precision = 1e-12, int max_iterations = 1000,
                   double deflation_tol = 1e-7) {
  assert(num_eigenvector >= 0);

  int iter = 1;
  auto mult = [&iter, &bonds, &couplings, &block](arma::Col<coeff_t> const &v,
                                                  arma::Col<coeff_t> &w) {
    auto ta = rightnow();
    Apply(bonds, couplings, block, v, block, w);
    Log(1, "Lanczos iteration {}", iter);
    timing(ta, rightnow(), "MVM", 1);
    ++iter;
  };

  auto converged = [num_eigenvector, precision](Tmatrix const &tmat) -> bool {
    return ConvergedEigenvalues(tmat, num_eigenvector, precision);
  };

  // First run for eigenvalues
  auto t0 = rightnow();
  // Allocate starting vector
  arma::Col<coeff_t> v0(block.size(), arma::fill::zeros);
  set_v0(v0);
  auto tmat = Lanczos(mult, v0, converged, max_iterations, deflation_tol);
  timing(t0, rightnow(), "Lanczos time (eigenvalue run)", 1);

  // Rerun for eigenvectors
  t0 = rightnow();
  iter = 1;
  auto tevecs = tmat.eigenvectors();
  arma::vec coefficients = tevecs.col(num_eigenvector);
  set_v0(v0);
  arma::Col<coeff_t> evec;
  std::tie(tmat, evec) = LanczosBuild(mult, v0, coefficients);
  timing(t0, rightnow(), "Lanczos time (eigenvector run)", 1);

  return {tmat, evec};
}

// Implementation with user-defined starting vector v0 (does a copy)
template <class coeff_t, class Block>
std::pair<Tmatrix, arma::Col<coeff_t>>
LanczosEigenvector(BondList const &bonds, Couplings const &couplings,
                   Block const &block, arma::Col<coeff_t> const &v0,
                   int num_eigenvector = 0, double precision = 1e-12,
                   int max_iterations = 1000, double deflation_tol = 1e-7) {

  auto set_v0 = [&v0](arma::Col<coeff_t> &v0_copy) { v0_copy = v0; };

  return LanczosEigenvector(bonds, couplings, block, set_v0, num_eigenvector,
                            precision, max_iterations, deflation_tol);
}

// Implementation with random real starting vector v0 (does NOT copy)
template <class Block>
std::pair<Tmatrix, arma::Col<double>>
LanczosEigenvectorReal(BondList const &bonds, Couplings const &couplings,
                       Block const &block, int num_eigenvector = 0,
                       double precision = 1e-12, int seed = 42,
                       int max_iterations = 1000, double deflation_tol = 1e-7) {

  // use different seeds for different MPI processes
  if constexpr (mpi::is_mpi_block<Block>) {
    seed += 0x01000193 * block.mpi_rank();
  }

  // Create random starting vector with normal distributed entries
  auto set_v0 = [&seed, &block](arma::Col<double> &v0) {
    uint32_t seed_modified = random::hash_combine(seed, random::hash(block));
    random::fill_random_normal_vector(v0, seed_modified);
  };

  // Run Lanczos algorithm
  return LanczosEigenvector<double, Block>(bonds, couplings, block, set_v0,
                                           num_eigenvector, precision,
                                           max_iterations, deflation_tol);
}

// Implementation with random complex starting vector v0 (does NOT copy)
template <class Block>
std::pair<Tmatrix, arma::Col<complex>>
LanczosEigenvectorCplx(BondList const &bonds, Couplings const &couplings,
                       Block const &block, int num_eigenvector = 0,
                       double precision = 1e-12, int seed = 42,
                       int max_iterations = 1000, double deflation_tol = 1e-7) {

  // use different seeds for different MPI processes
  if constexpr (mpi::is_mpi_block<Block>) {
    seed += 0x01000193 * block.mpi_rank();
  }

  // Create random starting vector with normal distributed entries
  auto set_v0 = [&seed, &block](arma::Col<complex> &v0) {
    uint32_t seed_modified = random::hash_combine(seed, random::hash(block));
    random::fill_random_normal_vector(v0, seed_modified);
  };

  // Run Lanczos algorithm
  return LanczosEigenvector<complex, Block>(bonds, couplings, block, set_v0,
                                            num_eigenvector, precision,
                                            max_iterations, deflation_tol);
}

} // namespace hydra
