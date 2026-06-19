// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "coeff.hpp"

#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>

#ifndef XDIAG_DISABLE_COLOR
#include <extern/fmt/color.hpp>
#endif

namespace xdiag {

Coeff::Coeff(std::string value) : value_(std::move(value)) {}
Coeff::Coeff(const char *value) : value_(std::string(value)) {}
Coeff::Coeff(double value) : value_(Scalar(value)) {}
Coeff::Coeff(complex value) : value_(Scalar(value)) {}
Coeff::Coeff(Scalar value) : value_(std::move(value)) {}

bool Coeff::isscalar() const { return std::holds_alternative<Scalar>(value_); }
bool Coeff::isstring() const {
  return std::holds_alternative<std::string>(value_);
}

Scalar Coeff::scalar() const try {
  if (const Scalar *v = std::get_if<Scalar>(&value_)) {
    return *v;
  } else {
    XDIAG_THROW("Cannot convert Coeff holding a string value to Scalar. "
                "Call OpSum::plain() first to resolve named coefficients.");
  }
}
XDIAG_CATCH

std::string Coeff::string() const try {
  if (const std::string *v = std::get_if<std::string>(&value_)) {
    return *v;
  } else {
    XDIAG_THROW("Cannot convert Coeff holding a Scalar value to std::string.");
  }
}
XDIAG_CATCH

bool Coeff::operator==(Coeff const &rhs) const { return value_ == rhs.value_; }
bool Coeff::operator!=(Coeff const &rhs) const { return !operator==(rhs); }

bool isscalar(Coeff const &c) { return c.isscalar(); }
bool isstring(Coeff const &c) { return c.isstring(); }
Scalar scalar(Coeff const &c) { return c.scalar(); }
std::string string(Coeff const &c) { return c.string(); }

Coeff operator*(Coeff const &lhs, Coeff const &rhs) try {
  if (lhs.isscalar() && rhs.isscalar()) {
    return Coeff(lhs.scalar() * rhs.scalar());
  } else {
    XDIAG_THROW("Cannot multiply Coefficients that contain string values. "
                "Call OpSum::plain() first to resolve named coefficients.");
  }
}
XDIAG_CATCH

Coeff operator*(Coeff const &lhs, Scalar const &rhs) try {
  if (lhs.isscalar()) {
    return Coeff(lhs.scalar() * rhs);
  } else {
    XDIAG_THROW("Cannot multiply a string Coeff with a Scalar. "
                "Call OpSum::plain() first to resolve named coefficients.");
  }
}
XDIAG_CATCH

Coeff operator*(Scalar const &lhs, Coeff const &rhs) try { return rhs * lhs; }
XDIAG_CATCH

std::ostream &operator<<(std::ostream &out, Coeff const &c) {
  if (c.isscalar()) {
    out << c.scalar();
  } else {
#ifndef XDIAG_DISABLE_COLOR
    out << fmt::format(fg(fmt::rgb(0x98C379)), "'{}'", c.string());
#else
    out << "'" << c.string() << "'";
#endif
  }
  return out;
}

std::string to_string(Coeff const &c) { return to_string_generic(c); }

} // namespace xdiag
