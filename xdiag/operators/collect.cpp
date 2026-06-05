// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "collect.hpp"

#include <map>
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

  // Fuse size-1 "Matrix" monomials acting on the same (equally ordered) site
  // vector by summing their coefficient-weighted matrices. This is valid for
  // the "Matrix" type alone: a Matrix Op is, by definition, the linear operator
  // given by its matrix on those sites, independent of any algebra, so
  //   c1 * Matrix(S, M1) + c2 * Matrix(S, M2) = Matrix(S, c1*M1 + c2*M2).
  std::map<SiteVec, Matrix> matrix_sums;
  OpSum result;
  for (auto const &[mono, c] : sums) {
    if (mono.size() == 1 && mono[0].type() == "Matrix" && mono[0].hassites() &&
        mono[0].hasmatrix()) {
      SiteVec const &sites = mono[0].sites();
      Matrix weighted = mono[0].matrix() * c;
      std::map<SiteVec, Matrix>::iterator it = matrix_sums.find(sites);
      if (it == matrix_sums.end()) {
        matrix_sums.emplace(sites, std::move(weighted));
      } else {
        it->second += weighted;
      }
    } else if (abs(c) > tol) {
      result += c * mono;
    }
  }

  // Emit the fused Matrix Ops, dropping any that summed to the zero matrix.
  for (auto const &[sites, m] : matrix_sums) {
    Matrix zero = arma::mat(m.n_rows(), m.n_cols(), arma::fill::zeros);
    if (!m.isapprox(zero, tol, tol)) {
      result += OpSum(Op("Matrix", sites, m));
    }
  }
  return result;
}
XDIAG_CATCH

} // namespace xdiag
