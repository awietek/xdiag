#pragma once

#include <tuple>
#include <utility>
#include <vector>

#include "extern/armadillo/armadillo"

#include <hydra/linalg/lanczos/lanczos_eigenvalues.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/states/state.h>
#include <hydra/utils/logger.h>

namespace hydra {

template <class Block>
double E0Real(BondList const &bondlist, Couplings const &couplings,
              Block &&block, double precision = 1e-12, int seed = 42,
              int max_iterations = 1000) {

  auto tmat = LanczosEigenvaluesReal(bondlist, couplings, block, 1, precision,
                                     seed, max_iterations);
  auto eigs = tmat.eigenvalues();
  if (eigs.size() == 0) {
    Log.err("Error: Tmatrix zero dimensional in E0Real");
    return 0.;
  } else {
    return tmat.eigenvalues()(0);
  }
}

template <class Block>
double E0Cplx(BondList const &bondlist, Couplings const &couplings,
              Block &&block, double precision = 1e-12, int seed = 42,
              int max_iterations = 1000) {
  auto tmat = LanczosEigenvaluesCplx(bondlist, couplings, block, 1, precision,
                                     seed, max_iterations);
  auto eigs = tmat.eigenvalues();
  if (eigs.size() == 0) {
    Log.err("Error: Tmatrix zero dimensional in E0Cplx");
    return 0.;
  } else {
    return tmat.eigenvalues()(0);
  }
}

template <class Block>
std::pair<double, StateReal<Block>>
GroundstateReal(BondList const &bondlist, Couplings const &couplings,
                Block const &block, double precision = 1e-12, int seed = 42,
                int max_iterations = 1000) {

  auto [tmat, v0] = LanczosEigenvectorReal(bondlist, couplings, block, 0,
                                           precision, seed, max_iterations);
  auto eigs = tmat.eigenvalues();
  if (eigs.size() == 0) {
    Log.err("Error: Tmatrix zero dimensional in GroundstateReal");
    return {0., State<double, Block>()};
  } else {
    return {tmat.eigenvalues()(0), State(block, v0)};
  }
}

template <class Block>
std::pair<double, StateCplx<Block>>
GroundstateCplx(BondList const &bondlist, Couplings const &couplings,
                Block const &block, double precision = 1e-12, int seed = 42,
                int max_iterations = 1000) {
  auto [tmat, v0] = LanczosEigenvectorCplx(bondlist, couplings, block, 0,
                                           precision, seed, max_iterations);
  auto eigs = tmat.eigenvalues();
  if (eigs.size() == 0) {
    Log.err("Error: Tmatrix zero dimensional in GroundstateCplx");
    return {0., State<complex, Block>()};
  } else {
    return {tmat.eigenvalues()(0), State(block, v0)};
  }
}

} // namespace hydra
