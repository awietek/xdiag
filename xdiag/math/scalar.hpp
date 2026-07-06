// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <variant>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Scalar is a type-erased numeric value: either double (real) or complex.
// It stays real until complex arithmetic forces a promotion; the variant
// then widens to complex and never narrows back.
//
// Use is<double>() / is<complex>() to query the active type.
// Use as<double>() to extract (throws if complex); as<complex>() always works.
// Use to_real() to narrow back to double if imag() is zero (throws otherwise).
class Scalar {
public:
  using value_t = std::variant<double, complex>;

  Scalar() = default;
  Scalar(double value);
  Scalar(complex value);

  Scalar(int64_t value); // prevents ambiguous integer literal overload
  Scalar(int value);

  bool operator==(Scalar const &rhs) const;
  bool operator!=(Scalar const &rhs) const;

  // Field operations. Mixed real/complex widens the result to complex.
  Scalar &operator+=(Scalar const &rhs);
  Scalar &operator-=(Scalar const &rhs);
  Scalar &operator*=(Scalar const &rhs);
  Scalar &operator/=(Scalar const &rhs);

  Scalar operator-() const;
  Scalar operator+(Scalar const &b) const;
  Scalar operator-(Scalar const &b) const;
  Scalar operator*(Scalar const &b) const;
  Scalar operator/(Scalar const &b) const;

  template <typename T> bool is() const; // is<double>() or is<complex>()
  template <typename T> T as() const;    // as<double>() throws if complex

  bool isreal() const;
  double real() const;
  double imag() const; // 0 for real Scalars
  double abs() const;
  Scalar conj() const; // identity for real Scalars
  Scalar to_real(
      double tol = 1e-12) const; // narrows to double; throws if |imag| > tol
  bool isapprox(Scalar const &y, double rtol, double atol) const;

private:
  value_t value_;
};

// Returns a zero of the same type (double or complex) as s.
Scalar zero(Scalar s);

bool isreal(Scalar const &s);
double real(Scalar const &s);
double imag(Scalar const &s);
double abs(Scalar const &s);
Scalar conj(Scalar const &s);
bool isapprox(Scalar const &a, Scalar const &b, double rtol = 1e-12,
              double atol = 1e-12);

std::ostream &operator<<(std::ostream &out, Scalar const &v);
std::string to_string(Scalar const &v);

} // namespace xdiag
