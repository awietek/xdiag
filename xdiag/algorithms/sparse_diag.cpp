// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "sparse_diag.hpp"

#include <xdiag/algorithms/lanczos/eigs_lanczos.hpp>
#include <xdiag/algorithms/lanczos/eigvals_lanczos.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag {

double eigval0(OpSum const &ops, Block const &block, double precision,
               int64_t max_iterations, int64_t random_seed) try {
  if (dim(block) == 0) {
    Log.warn("Warning: block zero dimensional in eigval0");
    return std::nan("");
  }
  auto res = eigvals_lanczos(ops, block, 1, precision, max_iterations, 1e-7,
                             random_seed);

  if (res.eigenvalues.size() == 0) {
    Log.warn("Warning: Tmatrix zero dimensional in eig0_cplx");
    return std::nan("");
  } else {
    return res.eigenvalues(0);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return std::nan("");
}

std::tuple<double, State> eig0(OpSum const &ops, Block const &block,
                               double precision, int64_t max_iterations,
                               int64_t random_seed) try {
  if (dim(block) == 0) {
    Log.warn("Warning: block zero dimensional in eigval0");
    return {std::nan(""), State()};
  }
  auto res =
      eigs_lanczos(ops, block, 1, precision, max_iterations, 1e-7, random_seed);

  if (res.eigenvalues.size() == 0) {
    Log.warn("Warning: Tmatrix zero dimensional in eig0_cplx");
    return {std::nan(""), State()};
  } else {
    return {res.eigenvalues(0), res.eigenvectors};
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return {std::nan(""), State()};
}

} // namespace xdiag
