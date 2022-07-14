#pragma once

#include <algorithm>

#include <hydra/linalg/lanczos/lanczos_convergence.h>
#include <hydra/linalg/lanczos/lanczos_generic.h>
#include <hydra/linalg/lanczos/tmatrix.h>

#include <hydra/mpi/dot_mpi.h>
#include <hydra/mpi/timing_mpi.h>

#include <hydra/blocks/blocks.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <lila/all.h>

namespace hydra {

// Implementation with a function setting the starting vector
template <class coeff_t, class Block, class set_v0_f>
std::pair<Tmatrix, lila::Vector<coeff_t>>
LanczosEigenvector(BondList const &bonds, Couplings const &couplings,
                   Block const &block, set_v0_f set_v0, int num_eigenvector = 0,
                   double precision = 1e-12, int max_iterations = 1000,
                   double deflation_tol = 1e-7) {

  using namespace lila;

  assert(num_eigenvector >= 0);

  // Allocate starting vector
  auto v0 = lila::Zeros<coeff_t>(block.size());

  // MPI Lanczos
#ifdef HYDRA_ENABLE_MPI
  if constexpr (detail::is_mpi_block<Block>) {

    int iter = 0;
    auto mult = [&iter, &bonds, &couplings, &block](
                    lila::Vector<coeff_t> const &v, lila::Vector<coeff_t> &w) {
      auto ta = rightnow_mpi();
      Apply(bonds, couplings, block, v, block, w);
      timing_mpi(ta, rightnow_mpi(), "MVM", 1);
      ++iter;
    };

    auto dot_mpi = [](lila::Vector<coeff_t> const &v,
                      lila::Vector<coeff_t> const &w) -> coeff_t {
      return DotMPI(v, w);
    };

    auto converged = [num_eigenvector, precision](Tmatrix const &tmat) -> bool {
      return ConvergedEigenvalues(tmat, num_eigenvector, precision);
    };

    // First run for eigenvalues
    auto t0 = rightnow_mpi();
    set_v0(v0);
    auto [tmat, vectors] =
        LanczosGeneric(mult, v0, dot_mpi, converged, lila::Vector<coeff_t>(),
                       max_iterations, deflation_tol);
    timing_mpi(t0, rightnow_mpi(), "Lanczos time (eigenvalue run)", 1);

    // Rerun for eigenvectors
    t0 = rightnow_mpi();
    iter = 0;
    auto tevecs = tmat.eigenvectors();
    auto coefficients = tevecs.col(num_eigenvector);
    set_v0(v0);
    std::tie(tmat, vectors) =
        LanczosGeneric(mult, v0, dot_mpi, converged, coefficients,
                       max_iterations, deflation_tol);
    timing_mpi(t0, rightnow_mpi(), "Lanczos time (eigenvector run)", 1);

    return {tmat, vectors};
  }

  // Serial Lanczos
  else {
#endif
    int iter = 0;
    auto mult = [&iter, &bonds, &couplings, &block](
                    lila::Vector<coeff_t> const &v, lila::Vector<coeff_t> &w) {
      auto ta = rightnow();
      Apply(bonds, couplings, block, v, block, w);
      timing(ta, rightnow(), "MVM", 1);
      ++iter;
    };

    auto dot = [](lila::Vector<coeff_t> const &v,
                  lila::Vector<coeff_t> const &w) -> coeff_t {
      return lila::Dot(v, w);
    };

    auto converged = [num_eigenvector, precision](Tmatrix const &tmat) -> bool {
      return ConvergedEigenvalues(tmat, num_eigenvector, precision);
    };

    // First run for eigenvalues
    auto t0 = rightnow();
    set_v0(v0);
    auto [tmat, evec] =
        LanczosGeneric(mult, v0, dot, converged, lila::Vector<coeff_t>(),
                       max_iterations, deflation_tol);
    timing(t0, rightnow(), "Lanczos time (first run)", 1);

    // Rerun for eigenvectors
    t0 = rightnow();
    iter = 0;
    auto tevecs = tmat.eigenvectors();
    auto coefficients = tevecs.col(num_eigenvector);
    set_v0(v0);
    std::tie(tmat, evec) = LanczosGeneric(
        mult, v0, dot, converged, coefficients, max_iterations, deflation_tol);
    timing(t0, rightnow(), "Lanczos time (eigenvector run)", 1);
    return {tmat, evec};
#ifdef HYDRA_ENABLE_MPI
  }
#endif
}

// Implementation with user-defined starting vector v0 (does a copy)
template <class coeff_t, class Block>
std::pair<Tmatrix, lila::Vector<coeff_t>>
LanczosEigenvector(BondList const &bonds, Couplings const &couplings,
                   Block const &block, lila::Vector<coeff_t> const &v0,
                   int num_eigenvector = 0, double precision = 1e-12,
                   int max_iterations = 1000, double deflation_tol = 1e-7) {

  auto set_v0 = [&v0](lila::Vector<coeff_t> &v0_copy) { v0_copy = v0; };

  return LanczosEigenvector(bonds, couplings, block, set_v0, num_eigenvector,
                            precision, max_iterations, deflation_tol);
}

// Implementation with random real starting vector v0 (does NOT copy)
template <class Block>
std::pair<Tmatrix, lila::Vector<double>>
LanczosEigenvectorReal(BondList const &bonds, Couplings const &couplings,
                       Block const &block, int num_eigenvector = 0,
                       double precision = 1e-12, int seed = 42,
                       int max_iterations = 1000, double deflation_tol = 1e-7) {

  using namespace lila;

  // use different seeds for different MPI processes
  if constexpr (detail::is_mpi_block<Block>) {
    seed += 0x01000193 * block.mpi_rank();
  }

  // Create random starting vector with normal distributed entries
  auto set_v0 = [&seed](lila::Vector<double> &v0) {
    normal_dist_t<double> dist(0., 1.);
    normal_gen_t<double> gen(dist, seed);
    Random(v0, gen);
  };

  // Run Lanczos algorithm
  return LanczosEigenvector<double, Block>(bonds, couplings, block, set_v0,
                                           num_eigenvector, precision,
                                           max_iterations, deflation_tol);
}

// Implementation with random complex starting vector v0 (does NOT copy)
template <class Block>
std::pair<Tmatrix, lila::Vector<complex>>
LanczosEigenvectorCplx(BondList const &bonds, Couplings const &couplings,
                       Block const &block, int num_eigenvector = 0,
                       double precision = 1e-12, int seed = 42,
                       int max_iterations = 1000, double deflation_tol = 1e-7) {

  using namespace lila;

  // use different seeds for different MPI processes
  if constexpr (detail::is_mpi_block<Block>) {
    seed += 0x01000193 * block.mpi_rank();
  }

  // Create random starting vector with normal distributed entries
  auto set_v0 = [&seed](lila::Vector<complex> &v0) {
    normal_dist_t<complex> dist(0., 1.);
    normal_gen_t<complex> gen(dist, seed);
    Random(v0, gen);
  };

  // Run Lanczos algorithm
  return LanczosEigenvector<complex, Block>(bonds, couplings, block, set_v0,
                                            num_eigenvector, precision,
                                            max_iterations, deflation_tol);
}

} // namespace hydra
