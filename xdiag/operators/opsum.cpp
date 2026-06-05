// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "opsum.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>

namespace xdiag {

// --- Constructors ---

OpSum::OpSum(Op const &op) : terms_({{Coeff(1.0), Monomial(op)}}) {}
OpSum::OpSum(Monomial const &mono) : terms_({{Coeff(1.0), mono}}) {}
OpSum::OpSum(Scalar const &coeff, Op const &op)
    : terms_({{Coeff(coeff), Monomial(op)}}) {}
OpSum::OpSum(Scalar const &coeff, Monomial const &mono)
    : terms_({{Coeff(coeff), mono}}) {}
OpSum::OpSum(Coeff const &coeff, Op const &op)
    : terms_({{coeff, Monomial(op)}}) {}
OpSum::OpSum(Coeff const &coeff, Monomial const &mono)
    : terms_({{coeff, mono}}) {}
OpSum::OpSum(std::string const &coeff, Op const &op)
    : OpSum(Coeff(coeff), op) {}
OpSum::OpSum(double coeff, Op const &op) : OpSum(Coeff(coeff), op) {}
OpSum::OpSum(complex coeff, Op const &op) : OpSum(Coeff(coeff), op) {}
OpSum::OpSum(std::string const &coeff, Monomial const &mono)
    : OpSum(Coeff(coeff), mono) {}
OpSum::OpSum(double coeff, Monomial const &mono) : OpSum(Coeff(coeff), mono) {}
OpSum::OpSum(complex coeff, Monomial const &mono) : OpSum(Coeff(coeff), mono) {}

// --- Scalar scaling ---

OpSum &OpSum::operator*=(double scalar) try {
  return operator*=(Scalar(scalar));
}
XDIAG_CATCH

OpSum &OpSum::operator*=(complex scalar) try {
  return operator*=(Scalar(scalar));
}
XDIAG_CATCH

OpSum &OpSum::operator*=(Scalar const &scalar) try {
  for (auto &[c, m] : terms_) {
    if (c.isscalar()) {
      c = Coeff(c.scalar() * scalar);
    } else {
      auto it = params_.find(c.string());
      if (it != params_.end()) {
        c = Coeff(it->second * scalar);
      } else {
        XDIAG_THROW(fmt::format(
            "Cannot scale OpSum: coefficient \"{}\" has not been defined. "
            "Call plain() first to resolve named coefficients.",
            c.string()));
      }
    }
  }
  return *this;
}
XDIAG_CATCH

OpSum &OpSum::operator/=(double scalar) try {
  return operator/=(Scalar(scalar));
}
XDIAG_CATCH

OpSum &OpSum::operator/=(complex scalar) try {
  return operator/=(Scalar(scalar));
}
XDIAG_CATCH

OpSum &OpSum::operator/=(Scalar const &scalar) try {
  return operator*=(Scalar(1.0) / scalar);
}
XDIAG_CATCH

// --- Addition / Subtraction ---

void OpSum::merge_params(std::map<std::string, Scalar> const &other) try {
  for (auto const &[key, val] : other) {
    auto it = params_.find(key);
    if (it == params_.end()) {
      params_[key] = val;
    } else if (it->second != val) {
      XDIAG_THROW(fmt::format(
          "Conflicting values for named coefficient \"{}\": {} vs {}", key,
          to_string(it->second), to_string(val)));
    }
  }
}
XDIAG_CATCH

OpSum &OpSum::operator+=(OpSum const &ops) try {
  terms_.insert(terms_.end(), ops.terms_.begin(), ops.terms_.end());
  merge_params(ops.params_);
  return *this;
}
XDIAG_CATCH

OpSum &OpSum::operator+=(Op const &op) try { return operator+=(OpSum(op)); }
XDIAG_CATCH

OpSum &OpSum::operator+=(Monomial const &mono) try {
  return operator+=(OpSum(mono));
}
XDIAG_CATCH

OpSum OpSum::operator+(OpSum const &ops) const try {
  OpSum result = *this;
  result += ops;
  return result;
}
XDIAG_CATCH

OpSum OpSum::operator+(Op const &op) const try { return operator+(OpSum(op)); }
XDIAG_CATCH

OpSum OpSum::operator+(Monomial const &mono) const try {
  return operator+(OpSum(mono));
}
XDIAG_CATCH

OpSum &OpSum::operator-=(OpSum const &ops) try {
  OpSum neg = ops;
  neg *= Scalar(-1.0);
  return operator+=(neg);
}
XDIAG_CATCH

OpSum &OpSum::operator-=(Op const &op) try { return operator-=(OpSum(op)); }
XDIAG_CATCH

OpSum &OpSum::operator-=(Monomial const &mono) try {
  return operator-=(OpSum(mono));
}
XDIAG_CATCH

OpSum OpSum::operator-(OpSum const &ops) const try {
  OpSum result = *this;
  result -= ops;
  return result;
}
XDIAG_CATCH

OpSum OpSum::operator-(Op const &op) const try { return operator-(OpSum(op)); }
XDIAG_CATCH

OpSum OpSum::operator-(Monomial const &mono) const try {
  return operator-(OpSum(mono));
}
XDIAG_CATCH

OpSum OpSum::operator-() const try {
  OpSum result = *this;
  result *= Scalar(-1.0);
  return result;
}
XDIAG_CATCH

// --- Algebra product ---

OpSum OpSum::operator*(OpSum const &rhs) const try {
  OpSum result;
  for (auto const &[cl, ml] : terms_) {
    for (auto const &[cr, mr] : rhs.terms_) {
      result.terms_.push_back({cl * cr, ml * mr});
    }
  }
  result.merge_params(params_);
  result.merge_params(rhs.params_);
  return result;
}
XDIAG_CATCH

OpSum &OpSum::operator*=(OpSum const &rhs) try {
  *this = operator*(rhs);
  return *this;
}
XDIAG_CATCH

OpSum OpSum::operator*(Op const &rhs) const try {
  return operator*(OpSum(rhs));
}
XDIAG_CATCH

OpSum &OpSum::operator*=(Op const &rhs) try { return operator*=(OpSum(rhs)); }
XDIAG_CATCH

OpSum OpSum::operator*(Monomial const &rhs) const try {
  return operator*(OpSum(rhs));
}
XDIAG_CATCH

OpSum &OpSum::operator*=(Monomial const &rhs) try {
  return operator*=(OpSum(rhs));
}
XDIAG_CATCH

// --- Named parameters ---

Scalar &OpSum::operator[](std::string const &name) { return params_[name]; }
Scalar const &OpSum::operator[](std::string const &name) const try {
  return params_.at(name);
}
XDIAG_CATCH

OpSum OpSum::plain() const try {
  OpSum result;
  for (auto const &[coeff, mono] : terms_) {
    if (coeff.isscalar()) {
      result.terms_.push_back({coeff, mono});
    } else {
      auto it = params_.find(coeff.string());
      if (it != params_.end()) {
        result.terms_.push_back({Coeff(it->second), mono});
      } else {
        XDIAG_THROW(fmt::format(
            "Cannot make OpSum plain: coefficient \"{}\" has not been defined.",
            coeff.string()));
      }
    }
  }
  return result;
}
XDIAG_CATCH

// --- Access ---

std::vector<Term> const &OpSum::terms() const noexcept { return terms_; }
std::map<std::string, Scalar> const &OpSum::params() const noexcept {
  return params_;
}
int64_t OpSum::size() const noexcept { return (int64_t)terms_.size(); }
OpSum::iterator_t OpSum::begin() const noexcept { return terms_.begin(); }
OpSum::iterator_t OpSum::end() const noexcept { return terms_.end(); }

bool OpSum::operator==(OpSum const &rhs) const {
  return (terms_ == rhs.terms_) && (params_ == rhs.params_);
}
bool OpSum::operator!=(OpSum const &rhs) const { return !operator==(rhs); }

bool OpSum::isreal() const try {
  // A string coefficient has no numeric value on its own; resolve named
  // coefficients to their defined values first (plain() throws if any named
  // coefficient is undefined).
  for (Term const &t : plain().terms_) {
    if (!xdiag::isreal(t.coeff.scalar())) {
      return false;
    }
    if (!t.monomial.isreal()) {
      return false;
    }
  }
  return true;
}
XDIAG_CATCH

bool isreal(OpSum const &ops) { return ops.isreal(); }

bool OpSum::empty() const { return terms_.empty(); }
bool empty(OpSum const &ops) { return ops.empty(); }

// --- Free operators ---

// Coeff/scalar * Op
OpSum operator*(double coeff, Op const &op) { return OpSum(Coeff(coeff), op); }
OpSum operator*(complex coeff, Op const &op) { return OpSum(Coeff(coeff), op); }
OpSum operator*(std::string const &coeff, Op const &op) {
  return OpSum(Coeff(coeff), op);
}
OpSum operator*(Coeff const &coeff, Op const &op) { return OpSum(coeff, op); }
OpSum operator*(Scalar const &coeff, Op const &op) { return OpSum(coeff, op); }
OpSum operator*(Op const &op, double coeff) { return coeff * op; }
OpSum operator*(Op const &op, complex coeff) { return coeff * op; }
OpSum operator*(Op const &op, std::string const &coeff) { return coeff * op; }
OpSum operator*(Op const &op, Coeff const &coeff) { return coeff * op; }
OpSum operator*(Op const &op, Scalar const &coeff) { return coeff * op; }

// Coeff/scalar * Monomial
OpSum operator*(double coeff, Monomial const &mono) {
  return OpSum(Coeff(coeff), mono);
}
OpSum operator*(complex coeff, Monomial const &mono) {
  return OpSum(Coeff(coeff), mono);
}
OpSum operator*(std::string const &coeff, Monomial const &mono) {
  return OpSum(Coeff(coeff), mono);
}
OpSum operator*(Coeff const &coeff, Monomial const &mono) {
  return OpSum(coeff, mono);
}
OpSum operator*(Scalar const &coeff, Monomial const &mono) {
  return OpSum(coeff, mono);
}

OpSum operator*(Monomial const &mono, double coeff) { return coeff * mono; }
OpSum operator*(Monomial const &mono, complex coeff) { return coeff * mono; }
OpSum operator*(Monomial const &mono, std::string const &coeff) {
  return coeff * mono;
}
OpSum operator*(Monomial const &mono, Coeff const &coeff) {
  return coeff * mono;
}
OpSum operator*(Monomial const &mono, Scalar const &coeff) {
  return coeff * mono;
}

// scalar * OpSum (scalar scaling)
OpSum operator*(double scalar, OpSum const &ops) {
  return Scalar(scalar) * ops;
}
OpSum operator*(complex scalar, OpSum const &ops) {
  return Scalar(scalar) * ops;
}
OpSum operator*(Scalar const &scalar, OpSum const &ops) {
  OpSum result = ops;
  result *= scalar;
  return result;
}
OpSum operator*(OpSum const &ops, double scalar) { return scalar * ops; }
OpSum operator*(OpSum const &ops, complex scalar) { return scalar * ops; }
OpSum operator*(OpSum const &ops, Scalar const &scalar) { return scalar * ops; }

OpSum operator/(OpSum const &ops, double scalar) {
  return ops / Scalar(scalar);
}
OpSum operator/(OpSum const &ops, complex scalar) {
  return ops / Scalar(scalar);
}
OpSum operator/(OpSum const &ops, Scalar const &scalar) {
  return ops * (Scalar(1.0) / scalar);
}

// Additive ops: Monomial on the left
OpSum operator+(Monomial const &lhs, Monomial const &rhs) {
  return OpSum(lhs) + OpSum(rhs);
}
OpSum operator-(Monomial const &lhs, Monomial const &rhs) {
  return OpSum(lhs) - OpSum(rhs);
}
OpSum operator-(Monomial const &mono) { return -OpSum(mono); }
OpSum operator+(Monomial const &lhs, OpSum const &rhs) {
  return OpSum(lhs) + rhs;
}
OpSum operator-(Monomial const &lhs, OpSum const &rhs) {
  return OpSum(lhs) - rhs;
}

// Algebra products: Op * OpSum, Monomial * OpSum
OpSum operator*(Op const &lhs, OpSum const &rhs) { return OpSum(lhs) * rhs; }
OpSum operator*(Monomial const &lhs, OpSum const &rhs) {
  return OpSum(lhs) * rhs;
}

// --- I/O ---

// Returns the visible (uncolored) width of a coefficient for alignment
// purposes. ANSI escape codes in the styled string are invisible, so we must
// not use the styled string's byte length for padding.
static size_t coeff_visible_width(Coeff const &c) {
  if (c.isscalar()) {
    return to_string(c.scalar()).size(); // scalar output has no ANSI codes
  } else {
    return c.string().size() + 2; // 'name' — name plus two single quotes
  }
}

std::ostream &operator<<(std::ostream &out, OpSum const &ops) {
  if (ops.size() == 0) {
    out << "0\n";
    return out;
  }

  // Collect styled coefficient strings, their visible widths, and mono strings
  std::vector<std::string> cstrs, mstrs;
  std::vector<size_t> cwidths;
  cstrs.reserve(ops.size());
  mstrs.reserve(ops.size());
  cwidths.reserve(ops.size());
  for (auto const &[coeff, mono] : ops) {
    cstrs.push_back(to_string(coeff));
    cwidths.push_back(coeff_visible_width(coeff));
    mstrs.push_back(to_string(mono));
  }

  // Max visible width determines the padding column
  size_t cw = *std::max_element(cwidths.begin(), cwidths.end());

  // Right-align using explicit space padding (not fmt width, which counts
  // bytes)
  for (size_t i = 0; i < cstrs.size(); ++i) {
    out << "  " << std::string(cw - cwidths[i], ' ') << cstrs[i] << "  "
        << mstrs[i] << "\n";
  }

  // Named parameters section
  if (!ops.params().empty()) {
    out << "\n";
    size_t nw = 0;
    for (auto const &[name, val] : ops.params())
      nw = std::max(nw, name.size());
    for (auto const &[name, val] : ops.params())
      out << fmt::format("  {0:>{1}} = {2}\n", name, nw, to_string(val));
  }
  return out;
}

std::string to_string(OpSum const &ops) { return to_string_generic(ops); }

} // namespace xdiag
