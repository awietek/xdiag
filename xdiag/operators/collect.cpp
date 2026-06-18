// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "collect.hpp"

#include <map>
#include <utility>
#include <vector>

#include <xdiag/math/matrix.hpp>
#include <xdiag/math/scalar.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

OpSum collect(OpSum const &ops, double tol) try {
  using SiteVec = std::vector<int64_t>;

  // Group terms by monomial, summing scalar coefficients.
  // std::map uses Monomial::operator< so the output is in canonical order.
  std::map<Monomial, Scalar> sums;
  for (auto const &[coeff, mono] : ops.plain()) {
    sums[mono] += coeff.scalar();
  }

  // Bucket size-1 "Matrix" monomials by their (equally ordered) site vector.
  // Several Matrix ops on the SAME sites get fused into one matrix:
  //   c1 * Matrix(S, M1) + c2 * Matrix(S, M2) = Matrix(S, c1*M1 + c2*M2),
  // which is valid because a Matrix Op is by definition the linear operator
  // given by its matrix, independent of any algebra. A LONE Matrix op on a site
  // vector is kept as c * Matrix(S, M) WITHOUT folding c into the matrix, so a
  // scalar coefficient (e.g. a momentum-symmetrization phase) stays visible to
  // callers that inspect coefficients, such as representation().
  std::map<SiteVec, std::vector<std::pair<Scalar, Matrix>>> matrix_buckets;
  OpSum result;
  for (auto const &[mono, c] : sums) {
    if (mono.size() == 1 && mono[0].type() == "Matrix" && mono[0].hassites() &&
        mono[0].hasmatrix()) {
      matrix_buckets[mono[0].sites()].push_back({c, mono[0].matrix()});
    } else if (abs(c) > tol) {
      result += c * mono;
    }
  }

  for (auto const &[sites, items] : matrix_buckets) {
    if (items.size() == 1) {
      Scalar const &c = items[0].first;
      if (abs(c) > tol) {
        result += c * OpSum(Op("Matrix", sites, items[0].second));
      }
    } else { // fuse: sum the coefficient-weighted matrices into one Matrix op
      Matrix fused = items[0].second * items[0].first;
      for (std::size_t i = 1; i < items.size(); ++i) {
        fused += items[i].second * items[i].first;
      }
      Matrix zero = arma::mat(fused.n_rows(), fused.n_cols(), arma::fill::zeros);
      if (!fused.isapprox(zero, tol, tol)) {
        result += OpSum(Op("Matrix", sites, fused));
      }
    }
  }
  return result;
}
XDIAG_CATCH

} // namespace xdiag
