#pragma once

#include <lila/all.h>
#include <tuple>
#include <vector>

#include <hydra/linalg/lanczos/lanczos_eigenvalues.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class Block>
double E0Real(BondList const &bondlist, Couplings const &couplings,
              Block &&block, double precision = 1e-12, int seed = 42,
              int max_iterations = 1000) {

  auto tmat = LanczosEigenvaluesReal(bondlist, couplings, block, 1, precision,
                                     seed, max_iterations);
  auto eigs = tmat.eigenvalues();
  if (eigs.size() == 0) {
    lila::Log.err("Tmatrix zero dimensional in E0Real");
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
    lila::Log.err("Tmatrix zero dimensional in E0Cplx");
    return 0.;
  } else {
    return tmat.eigenvalues()(0);
  }
}

} // namespace hydra
