// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "monomial.hpp"

#include <xdiag/operators/hc.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/to_string_generic.hpp>

namespace xdiag {

Monomial::Monomial(Op const &op) : ops_({op}) {}
Monomial::Monomial(std::initializer_list<Op> ops) : ops_(ops) {}
Monomial::Monomial(std::vector<Op> const &ops) : ops_(ops) {}

int64_t Monomial::size() const noexcept { return (int64_t)ops_.size(); }
bool Monomial::empty() const noexcept { return ops_.empty(); }

Op const &Monomial::operator[](int64_t idx) const try { return ops_.at(idx); }
XDIAG_CATCH

std::vector<Op> const &Monomial::ops() const noexcept { return ops_; }

Monomial::iterator_t Monomial::begin() const noexcept { return ops_.begin(); }
Monomial::iterator_t Monomial::end() const noexcept { return ops_.end(); }

Monomial Monomial::operator*(Op const &rhs) const {
  Monomial result = *this;
  result.ops_.push_back(rhs);
  return result;
}

Monomial Monomial::operator*(Monomial const &rhs) const {
  Monomial result = *this;
  result.ops_.insert(result.ops_.end(), rhs.ops_.begin(), rhs.ops_.end());
  return result;
}

Monomial &Monomial::operator*=(Op const &rhs) {
  ops_.push_back(rhs);
  return *this;
}

Monomial &Monomial::operator*=(Monomial const &rhs) {
  ops_.insert(ops_.end(), rhs.ops_.begin(), rhs.ops_.end());
  return *this;
}

bool Monomial::operator==(Monomial const &rhs) const noexcept {
  return ops_ == rhs.ops_;
}
bool Monomial::operator!=(Monomial const &rhs) const noexcept {
  return !operator==(rhs);
}
bool Monomial::operator<(Monomial const &rhs) const noexcept {
  // Shorter monomials come first; equal length: lexicographic by Op::operator<
  if (ops_.size() != rhs.ops_.size())
    return ops_.size() < rhs.ops_.size();
  for (int64_t i = 0; i < (int64_t)ops_.size(); ++i) {
    if (ops_[i] != rhs.ops_[i])
      return ops_[i] < rhs.ops_[i];
  }
  return false; // equal
}

Monomial Monomial::hc() const {
  Monomial result;
  result.ops_.reserve(ops_.size());
  for (int64_t i = (int64_t)ops_.size() - 1; i >= 0; --i) {
    result.ops_.push_back(xdiag::hc(ops_[i]));
  }
  return result;
}

bool Monomial::isreal() const {
  for (Op const &op : ops_) {
    if (!op.isreal()) {
      return false;
    }
  }
  return true;
}

bool isreal(Monomial const &m) { return m.isreal(); }

Monomial operator*(Op const &lhs, Op const &rhs) {
  return Monomial({lhs, rhs});
}

Monomial operator*(Op const &lhs, Monomial const &rhs) {
  Monomial result(lhs);
  result *= rhs;
  return result;
}

std::ostream &operator<<(std::ostream &out, Monomial const &m) {
  if (m.empty()) {
    out << "1";
  } else {
    for (int64_t i = 0; i < m.size(); ++i) {
      if (i > 0)
        out << " * ";
      out << m[i];
    }
  }
  return out;
}

std::string to_string(Monomial const &m) { return to_string_generic(m); }

bool Term::operator==(Term const &rhs) const noexcept {
  return coeff == rhs.coeff && monomial == rhs.monomial;
}
bool Term::operator!=(Term const &rhs) const noexcept {
  return !operator==(rhs);
}

} // namespace xdiag
