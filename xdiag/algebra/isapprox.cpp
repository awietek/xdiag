// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "isapprox.hpp"

#include <cmath>
#include <vector>

#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/math/matrix.hpp>
#include <xdiag/math/scalar.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

bool isapprox(Op const &op1, Op const &op2, double rtol, double atol) {
  if (op1.type() != op2.type()) {
    return false;
  } else if (op1.hassites() != op2.hassites()) {
    return false;
  } else if (op1.hassites() && op1.sites() != op2.sites()) {
    return false;
  } else if (op1.hasmatrix() != op2.hasmatrix()) {
    return false;
  } else if (op1.hasmatrix()) {
    return xdiag::isapprox(op1.matrix(), op2.matrix(), rtol, atol);
  }
  return true;
}

namespace {

// A term of a normal-ordered OpSum: its scalar coefficient and its monomial.
using OpTerm = std::pair<Scalar, Monomial>;

// Two monomials have the same "skeleton" if they are the same product of
// operators ignoring any Matrix VALUES: same length, and op-by-op the same type
// and sites. We never order or compare matrices here -- a Matrix has no
// sensible (and no danger-free) operator<, and the (coefficient, matrix) split
// is not unique since c*Matrix(M) == (c*a)*Matrix(M/a). Whether two
// same-skeleton terms agree (up to a scalar) is decided afterwards by isapprox
// on the coefficient-weighted matrices.
bool same_skeleton(Monomial const &a, Monomial const &b) {
  if (a.size() != b.size()) {
    return false;
  }
  for (int64_t i = 0; i < a.size(); ++i) {
    if (a[i].type() != b[i].type()) {
      return false;
    }
    if (a[i].hassites() != b[i].hassites()) {
      return false;
    }
    if (a[i].hassites() && a[i].sites() != b[i].sites()) {
      return false;
    }
  }
  return true;
}

bool is_matrix_mono(Monomial const &m) {
  return m.size() == 1 && m[0].type() == "Matrix" && m[0].hasmatrix();
}

// Is the coefficient-weighted operator of this term (approximately) zero?
bool is_zero(OpTerm const &t, double rtol, double atol) {
  if (is_matrix_mono(t.second)) {
    Matrix v = t.second[0].matrix() * t.first;
    Matrix zero = arma::mat(v.n_rows(), v.n_cols(), arma::fill::zeros);
    return v.isapprox(zero, rtol, atol);
  }
  return std::abs(t.first.as<complex>()) <= atol;
}

// Does term `a` equal `scale * b` ? Requires the same skeleton, then compares
// the coefficient-weighted matrices (resp. coefficients) with isapprox.
bool term_isapprox_scaled(OpTerm const &a, OpTerm const &b, complex scale,
                          double rtol, double atol) {
  if (!same_skeleton(a.second, b.second)) {
    return false;
  }
  if (is_matrix_mono(a.second)) {
    Matrix av = a.second[0].matrix() * a.first;
    Matrix bv = b.second[0].matrix() * (b.first * Scalar(scale));
    return xdiag::isapprox(av, bv, rtol, atol);
  }
  complex ca = a.first.as<complex>();
  complex cb = scale * b.first.as<complex>();
  return std::abs(ca - cb) <= atol + rtol * std::abs(cb);
}

// The scalar lambda with (coeff-weighted operator of a) = lambda * (that of b),
// read off the dominant component. Assumes same_skeleton(a, b) and b nonzero.
complex term_ratio(OpTerm const &a, OpTerm const &b) {
  if (is_matrix_mono(a.second)) {
    Matrix ma = a.second[0].matrix() * a.first;
    Matrix mb = b.second[0].matrix() * b.first;
    arma::cx_vec va =
        ma.isreal() ? arma::conv_to<arma::cx_vec>::from(
                          arma::vectorise(ma.as<arma::mat>()))
                    : arma::cx_vec(arma::vectorise(ma.as<arma::cx_mat>()));
    arma::cx_vec vb =
        mb.isreal() ? arma::conv_to<arma::cx_vec>::from(
                          arma::vectorise(mb.as<arma::mat>()))
                    : arma::cx_vec(arma::vectorise(mb.as<arma::cx_mat>()));
    arma::uword idx = arma::abs(vb).index_max();
    return va[idx] / vb[idx];
  }
  return a.first.as<complex>() / b.first.as<complex>();
}

// Non-zero terms of an OpSum brought into normal order.
std::vector<OpTerm> nonzero_terms(OpSum const &ops, algebra::Algebra const &alg,
                                  double rtol, double atol) {
  std::vector<OpTerm> terms;
  for (auto const &[c, mono] : algebra::normal_order(ops, alg).plain()) {
    OpTerm t{c.scalar(), mono};
    if (!is_zero(t, rtol, atol)) {
      terms.push_back(t);
    }
  }
  return terms;
}

// Does A equal `scale * B` term-by-term? Linear search: each A term must pair
// with a distinct B term via term_isapprox_scaled.
bool lists_match_scaled(std::vector<OpTerm> const &A,
                        std::vector<OpTerm> const &B, complex scale,
                        double rtol, double atol) {
  if (A.size() != B.size()) {
    return false;
  }
  std::vector<bool> used(B.size(), false);
  for (OpTerm const &a : A) {
    bool found = false;
    for (std::size_t j = 0; j < B.size(); ++j) {
      if (!used[j] && term_isapprox_scaled(a, B[j], scale, rtol, atol)) {
        used[j] = true;
        found = true;
        break;
      }
    }
    if (!found) {
      return false;
    }
  }
  return true;
}

} // namespace

bool isapprox(OpSum const &ops1, OpSum const &ops2,
              algebra::Algebra const &algebra, double rtol, double atol) try {
  std::vector<OpTerm> t1 = nonzero_terms(ops1, algebra, rtol, atol);
  std::vector<OpTerm> t2 = nonzero_terms(ops2, algebra, rtol, atol);
  return lists_match_scaled(t1, t2, complex(1.0, 0.0), rtol, atol);
}
XDIAG_CATCH

std::optional<Scalar> isapprox_multiple(OpSum const &ops1, OpSum const &ops2,
                                        algebra::Algebra const &algebra,
                                        double rtol, double atol) try {
  std::vector<OpTerm> t1 = nonzero_terms(ops1, algebra, rtol, atol);
  std::vector<OpTerm> t2 = nonzero_terms(ops2, algebra, rtol, atol);

  if (t1.empty() && t2.empty()) {
    return Scalar(0.0);
  }
  if (t1.empty() || t2.empty() || (t1.size() != t2.size())) {
    return std::nullopt;
  }

  // Candidate lambda: pair the first term of ops2 with the same-skeleton term
  // of ops1, then read off the ratio of their coefficient-weighted operators.
  OpTerm const &ref = t2.front();
  OpTerm const *match = nullptr;
  for (OpTerm const &a : t1) {
    if (same_skeleton(a.second, ref.second)) {
      match = &a;
      break;
    }
  }
  if (match == nullptr) {
    return std::nullopt;
  }
  complex lambda = term_ratio(*match, ref);

  // Verify ops1 = lambda * ops2 on every term.
  if (!lists_match_scaled(t1, t2, lambda, rtol, atol)) {
    return std::nullopt;
  }
  // Return a real Scalar when lambda is real within tolerance, so that callers
  // that expect a real result (e.g. a real correlation matrix) and callers that
  // inspect Scalar::isreal() (e.g. representation() deciding real vs complex
  // characters) behave as they did with exact arithmetic.
  if (std::abs(lambda.imag()) <= atol + rtol * std::abs(lambda)) {
    return Scalar(lambda.real());
  }
  return Scalar(lambda);
}
XDIAG_CATCH

} // namespace xdiag
