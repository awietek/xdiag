// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <variant>

#include <xdiag/math/complex.hpp>
#include <xdiag/math/scalar.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Coeff represents an operator coefficient that is either a numeric Scalar
// (double or complex) or a named string constant to be resolved later via
// OpSum::plain().
class Coeff {
public:
  using value_t = std::variant<Scalar, std::string>;

  Coeff() = default;
  explicit Coeff(std::string value);
  explicit Coeff(const char *value);
  explicit Coeff(double value);
  explicit Coeff(complex value);
  explicit Coeff(Scalar value);

  bool operator==(Coeff const &rhs) const;
  bool operator!=(Coeff const &rhs) const;

  bool isscalar() const;
  bool isstring() const;
  Scalar scalar() const;
  std::string string() const;

private:
  value_t value_;
};

bool isscalar(Coeff const &c);
bool isstring(Coeff const &c);
Scalar scalar(Coeff const &c);
std::string string(Coeff const &c);

// Multiply two Coeffs. Both must be scalar (throws if either is a string).
// Call OpSum::plain() first to resolve named coefficients before multiplying.
Coeff operator*(Coeff const &lhs, Coeff const &rhs);
Coeff operator*(Coeff const &lhs, Scalar const &rhs);
Coeff operator*(Scalar const &lhs, Coeff const &rhs);

std::ostream &operator<<(std::ostream &out, Coeff const &c);
std::string to_string(Coeff const &c);

} // namespace xdiag
