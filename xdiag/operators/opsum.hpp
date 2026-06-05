// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <map>
#include <string>
#include <vector>

#include <xdiag/math/scalar.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// OpSum represents a linear combination of Monomials (ordered products of Ops)
// with Coeff coefficients. This is the algebra structure:
//   OpSum = sum_i { coeff_i * (A_{i,1} * A_{i,2} * ... * A_{i,n_i}) }
//
// Vector space operations (+, -, scalar*) combine terms additively.
// Algebra product (*) forms the non-commutative product distributively:
//   (sum_i c_i M_i) * (sum_j d_j N_j) = sum_{i,j} (c_i * d_j) (M_i ++ N_j)
//
// Named string coefficients can be resolved to scalars via plain().
class OpSum {
public:
  using iterator_t = std::vector<Term>::const_iterator;

  XDIAG_API OpSum() = default;
  XDIAG_API explicit OpSum(Op const &op);
  XDIAG_API explicit OpSum(Monomial const &mono);
  XDIAG_API OpSum(Scalar const &coeff, Op const &op);
  XDIAG_API OpSum(Scalar const &coeff, Monomial const &mono);
  XDIAG_API OpSum(Coeff const &coeff, Op const &op);
  XDIAG_API OpSum(Coeff const &coeff, Monomial const &mono);
  XDIAG_API OpSum(std::string const &coeff, Op const &op);
  XDIAG_API OpSum(double coeff, Op const &op);
  XDIAG_API OpSum(complex coeff, Op const &op);
  XDIAG_API OpSum(std::string const &coeff, Monomial const &mono);
  XDIAG_API OpSum(double coeff, Monomial const &mono);
  XDIAG_API OpSum(complex coeff, Monomial const &mono);

  // Vector space: addition and subtraction
  XDIAG_API OpSum &operator+=(OpSum const &ops);
  XDIAG_API OpSum &operator+=(Op const &op);
  XDIAG_API OpSum &operator+=(Monomial const &mono);
  XDIAG_API OpSum operator+(OpSum const &ops) const;
  XDIAG_API OpSum operator+(Op const &op) const;
  XDIAG_API OpSum operator+(Monomial const &mono) const;
  XDIAG_API OpSum &operator-=(OpSum const &ops);
  XDIAG_API OpSum &operator-=(Op const &op);
  XDIAG_API OpSum &operator-=(Monomial const &mono);
  XDIAG_API OpSum operator-(OpSum const &ops) const;
  XDIAG_API OpSum operator-(Op const &op) const;
  XDIAG_API OpSum operator-(Monomial const &mono) const;
  XDIAG_API OpSum operator-() const;

  // Scalar scaling (multiplies each coefficient by the scalar)
  XDIAG_API OpSum &operator*=(double scalar);
  XDIAG_API OpSum &operator*=(complex scalar);
  XDIAG_API OpSum &operator/=(double scalar);
  XDIAG_API OpSum &operator/=(complex scalar);
  OpSum &operator*=(Scalar const &scalar);
  OpSum &operator/=(Scalar const &scalar);

  // Algebra product (non-commutative, distributes over addition)
  XDIAG_API OpSum operator*(OpSum const &rhs) const;
  XDIAG_API OpSum &operator*=(OpSum const &rhs);
  XDIAG_API OpSum operator*(Op const &rhs) const;
  XDIAG_API OpSum &operator*=(Op const &rhs);
  XDIAG_API OpSum operator*(Monomial const &rhs) const;
  XDIAG_API OpSum &operator*=(Monomial const &rhs);

  // Named parameters: ops["J"] = 1.5; ops.plain() resolves them
  XDIAG_API Scalar &operator[](std::string const &name);
  XDIAG_API Scalar const &operator[](std::string const &name) const;
  XDIAG_API OpSum plain() const;

  // Access
  XDIAG_API std::vector<Term> const &terms() const noexcept;
  XDIAG_API std::map<std::string, Scalar> const &params() const noexcept;
  XDIAG_API int64_t size() const noexcept;
  XDIAG_API iterator_t begin() const noexcept;
  XDIAG_API iterator_t end() const noexcept;

  XDIAG_API bool operator==(OpSum const &rhs) const;
  XDIAG_API bool operator!=(OpSum const &rhs) const;

  XDIAG_API bool isreal() const;
  XDIAG_API bool empty() const;

private:
  std::vector<Term> terms_;
  std::map<std::string, Scalar> params_;

  void merge_params(std::map<std::string, Scalar> const &other);
};

// --- Free operators ---

// Coeff/scalar * Op, Op * Coeff/scalar  →  single-term OpSum
XDIAG_API OpSum operator*(double coeff, Op const &op);
XDIAG_API OpSum operator*(complex coeff, Op const &op);
XDIAG_API OpSum operator*(std::string const &coeff, Op const &op);
XDIAG_API OpSum operator*(Coeff const &coeff, Op const &op);
XDIAG_API OpSum operator*(Scalar const &coeff, Op const &op);
XDIAG_API OpSum operator*(Op const &op, double coeff);
XDIAG_API OpSum operator*(Op const &op, complex coeff);
XDIAG_API OpSum operator*(Op const &op, std::string const &coeff);
XDIAG_API OpSum operator*(Op const &op, Coeff const &coeff);
XDIAG_API OpSum operator*(Op const &op, Scalar const &coeff);

// Coeff/scalar * Monomial, Monomial * Coeff/scalar  →  single-term OpSum
XDIAG_API OpSum operator*(double coeff, Monomial const &mono);
XDIAG_API OpSum operator*(complex coeff, Monomial const &mono);
XDIAG_API OpSum operator*(std::string const &coeff, Monomial const &mono);
XDIAG_API OpSum operator*(Coeff const &coeff, Monomial const &mono);
XDIAG_API OpSum operator*(Scalar const &coeff, Monomial const &mono);
XDIAG_API OpSum operator*(Monomial const &mono, double coeff);
XDIAG_API OpSum operator*(Monomial const &mono, complex coeff);
XDIAG_API OpSum operator*(Monomial const &mono, std::string const &coeff);
XDIAG_API OpSum operator*(Monomial const &mono, Coeff const &coeff);
XDIAG_API OpSum operator*(Monomial const &mono, Scalar const &coeff);

// scalar * OpSum, OpSum * scalar  →  scaled copy
XDIAG_API OpSum operator*(double scalar, OpSum const &ops);
XDIAG_API OpSum operator*(complex scalar, OpSum const &ops);
XDIAG_API OpSum operator*(Scalar const &scalar, OpSum const &ops);
XDIAG_API OpSum operator*(OpSum const &ops, double scalar);
XDIAG_API OpSum operator*(OpSum const &ops, complex scalar);
XDIAG_API OpSum operator*(OpSum const &ops, Scalar const &scalar);
XDIAG_API OpSum operator/(OpSum const &ops, double scalar);
XDIAG_API OpSum operator/(OpSum const &ops, complex scalar);
XDIAG_API OpSum operator/(OpSum const &ops, Scalar const &scalar);

// Additive ops with Monomial on the left (Monomial → OpSum implicit conversion
// only works for member operator+, not when Monomial is the left operand)
XDIAG_API OpSum operator+(Monomial const &lhs, Monomial const &rhs);
XDIAG_API OpSum operator-(Monomial const &lhs, Monomial const &rhs);
XDIAG_API OpSum operator-(Monomial const &mono); // unary negation
XDIAG_API OpSum operator+(Monomial const &lhs, OpSum const &rhs);
XDIAG_API OpSum operator-(Monomial const &lhs, OpSum const &rhs);

// Algebra products: Op * OpSum, Monomial * OpSum
XDIAG_API OpSum operator*(Op const &lhs, OpSum const &rhs);
XDIAG_API OpSum operator*(Monomial const &lhs, OpSum const &rhs);

XDIAG_API bool isreal(OpSum const &ops);
XDIAG_API bool empty(OpSum const &ops);

XDIAG_API std::ostream &operator<<(std::ostream &out, OpSum const &ops);
XDIAG_API std::string to_string(OpSum const &ops);

} // namespace xdiag
