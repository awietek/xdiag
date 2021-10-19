#pragma once

#include <hydra/linalg/lanczos/lanczos_convergence.h>
#include <hydra/linalg/lanczos/lanczos_generic.h>
#include <hydra/linalg/lanczos/tmatrix.h>

#include <hydra/mpi/dot_mpi.h>
#include <hydra/mpi/timing_mpi.h>

#include <hydra/models/models.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <lila/all.h>

namespace hydra {

// Lanczos which  overwrites starting vector v0
template <class coeff_t, class Block>
Tmatrix LanczosEigenvaluesInplace(
    BondList const &bonds, Couplings const &couplings, Block const &block,
    lila::Vector<coeff_t> &v0, int num_eigenvalue = 0, double precision = 1e-12,
    int max_iterations = 1000, double deflation_tol = 1e-7) {

  using namespace lila;

  // MPI Lanczos
  if constexpr (is_mpi_block<Block>) {

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

    auto converged = [num_eigenvalue, precision](Tmatrix const &tmat) -> bool {
      return ConvergedEigenvalues(tmat, num_eigenvalue, precision);
    };

    auto t0 = rightnow_mpi();
    auto [tmat, vectors] =
        LanczosGeneric(mult, v0, dot_mpi, converged, lila::Matrix<coeff_t>(),
                       max_iterations, deflation_tol);
    (void)vectors;
    timing_mpi(t0, rightnow_mpi(), "Lanczos time", 1);
    return tmat;
  }

  // Serial Lanczos
  else {

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

    auto converged = [num_eigenvalue, precision](Tmatrix const &tmat) -> bool {
      return ConvergedEigenvalues(tmat, num_eigenvalue, precision);
    };

    auto t0 = rightnow();
    auto [tmat, vectors] =
        LanczosGeneric(mult, v0, dot, converged, lila::Matrix<coeff_t>(),
                       max_iterations, deflation_tol);
    (void)vectors;
    timing(t0, rightnow(), "Lanczos time", 1);
    return tmat;
  }
}

// Lanczos which does not overwrite v0
template <class coeff_t, class Block>
Tmatrix LanczosEigenvalues(BondList const &bonds, Couplings const &couplings,
                           Block const &block, lila::Vector<coeff_t> v0,
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

  using namespace lila;

  // use different seeds for different MPI processes
  if constexpr (is_mpi_block<Block>) {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    seed += 0x01000193 * mpi_rank;
  }

  // Create random starting vector with normal distributed entries
  normal_dist_t<double> dist(0., 1.);
  normal_gen_t<double> gen(dist, seed);
  auto v0 = Random(block.size(), gen);

  // Run Lanczos algorithm
  return LanczosEigenvaluesInplace(bonds, couplings, block, v0, num_eigenvalue,
                                   precision, max_iterations, deflation_tol);
}

template <class Block>
Tmatrix LanczosEigenvaluesCplx(BondList const &bonds,
                               Couplings const &couplings, Block const &block,
                               int num_eigenvalue = 0, double precision = 1e-12,
                               int seed = 42, int max_iterations = 1000,
                               double deflation_tol = 1e-7) {

  using namespace lila;

  // use different seeds for different MPI processes
  if constexpr (is_mpi_block<Block>) {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    seed += 0x01000193 * mpi_rank;
  }

  // Create random starting vector with normal distributed entries
  normal_dist_t<complex> dist(0., 1.);
  normal_gen_t<complex> gen(dist, seed);
  auto v0 = Random(block.size(), gen);

  // Run Lanczos algorithm
  return LanczosEigenvaluesInplace(bonds, couplings, block, v0, num_eigenvalue,
                                   precision, max_iterations, deflation_tol);
}

} // namespace hydra
