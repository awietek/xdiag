#pragma once

#include <tuple>
#include <utility>
#include <vector>

#include "extern/armadillo/armadillo"

#include <hydra/algorithms/lanczos/lanczos_eigenvalues.h>
#include <hydra/algorithms/lanczos/lanczos_eigenvector.h>
#include <hydra/operators/bondlist.h>
#include <hydra/states/state.h>
#include <hydra/utils/logger.h>

namespace hydra {

double eig0_real(BondList const &bondlist, Block const &block,
                 double precision = 1e-12, int seed = 42,
                 int max_iterations = 1000) {

  auto tmat = lanczos_eigenvalues_real(bondlist, block, 1, precision, seed,
                                       max_iterations);

  auto eigs = tmat.eigenvalues();
  if (eigs.size() == 0) {
    Log.err("Error: Tmatrix zero dimensional in eig0_real");
    return 0.;
  } else {
    return tmat.eigenvalues()(0);
  }
}

double eig0_cplx(BondList const &bondlist, Block const &block,
                 double precision = 1e-12, int seed = 42,
                 int max_iterations = 1000) {
  auto tmat = lanczos_eigenvalues_cplx(bondlist, block, 1, precision, seed,
                                       max_iterations);
  auto eigs = tmat.eigenvalues();
  if (eigs.size() == 0) {
    Log.err("Error: Tmatrix zero dimensional in eig0_cplx");
    return 0.;
  } else {
    return tmat.eigenvalues()(0);
  }
}

double eig0(BondList const &bondlist, Block const &block,
            double precision = 1e-12, int seed = 42,
            int max_iterations = 1000) {
  return eig0_cplx(bondlist, block, precision, seed, max_iterations);
}

StateReal groundstate_real(BondList const &bondlist, Block const &block,
                           double precision = 1e-12, int seed = 42,
                           int max_iterations = 1000) {

  auto [tmat, v0] = lanczos_eigenvector_real(bondlist, block, 0, precision,
                                             seed, max_iterations);
  auto eigs = tmat.eigenvalues();
  if (eigs.size() == 0) {
    Log.err("Error: Tmatrix zero dimensional in groundstate_real");
    return State<double>();
  } else {
    return State(block, v0);
  }
}

StateCplx groundstate_cplx(BondList const &bondlist, Block const &block,
                           double precision = 1e-12, int seed = 42,
                           int max_iterations = 1000) {
  auto [tmat, v0] = lanczos_eigenvector_cplx(bondlist, block, 0, precision,
                                             seed, max_iterations);
  auto eigs = tmat.eigenvalues();
  if (eigs.size() == 0) {
    Log.err("Error: Tmatrix zero dimensional in groundstate_cplx");
    return State<complex>();
  } else {
    return State(block, v0);
  }
}

StateCplx groundstate(BondList const &bondlist, Block const &block,
                      double precision = 1e-12, int seed = 42,
                      int max_iterations = 1000) {
  return groundstate_cplx(bondlist, block, precision, seed, max_iterations);
}

std::pair<double, StateReal> eig0_groundstate_real(BondList const &bondlist,
                                                   Block const &block,
                                                   double precision = 1e-12,
                                                   int seed = 42,
                                                   int max_iterations = 1000) {

  auto [tmat, v0] = lanczos_eigenvector_real(bondlist, block, 0, precision,
                                             seed, max_iterations);
  auto eigs = tmat.eigenvalues();
  if (eigs.size() == 0) {
    Log.err("Error: Tmatrix zero dimensional in groundstate_real");
    return {0., State<double>()};
  } else {
    return {tmat.eigenvalues()(0), State(block, v0)};
  }
}

std::pair<double, StateCplx> eig0_groundstate_cplx(BondList const &bondlist,
                                                   Block const &block,
                                                   double precision = 1e-12,
                                                   int seed = 42,
                                                   int max_iterations = 1000) {
  auto [tmat, v0] = lanczos_eigenvector_cplx(bondlist, block, 0, precision,
                                             seed, max_iterations);
  auto eigs = tmat.eigenvalues();
  if (eigs.size() == 0) {
    Log.err("Error: Tmatrix zero dimensional in groundstate_cplx");
    return {0., State<complex>()};
  } else {
    return {tmat.eigenvalues()(0), State(block, v0)};
  }
}

std::pair<double, StateCplx> eig0_groundstate(BondList const &bondlist,
                                              Block const &block,
                                              double precision = 1e-12,
                                              int seed = 42,
                                              int max_iterations = 1000) {
  return eig0_groundstate_cplx(bondlist, block, precision, seed,
                               max_iterations);
}

} // namespace hydra