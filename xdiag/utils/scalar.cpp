// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "scalar.hpp"

namespace xdiag {

Scalar::Scalar(double value) : value_(value) {}
Scalar::Scalar(complex value) : value_(value) {}

bool Scalar::isreal() const { return std::holds_alternative<double>(value_); }

double Scalar::real() const {
  return std::visit(overload{[](double const &a) { return a; },
                             [](complex const &a) { return std::real(a); }},
                    value_);
}

double Scalar::imag() const {
  return std::visit(overload{[](double const &a) { return 0.0; },
                             [](complex const &a) { return std::imag(a); }},
                    value_);
}

Scalar Scalar::conj() const {
  return std::visit(
      overload{[](double const &a) { return Scalar(a); },
               [](complex const &a) { return Scalar(std::conj(a)); }},
      value_);
}

double Scalar::abs() const {
  return std::visit([](auto &&a) { return std::abs(a); }, value_);
}

bool Scalar::isapprox(Scalar const &y, double rtol = 1e-12,
                      double atol = 1e-12) const {
  Scalar const &x = *this;
  return (x - y).abs() <= (atol + rtol * y.abs());
}

bool Scalar::operator==(Scalar const &rhs) const {
  return value_ == rhs.value_;
}

bool Scalar::operator!=(Scalar const &rhs) const { return !operator==(rhs); }

Scalar &Scalar::operator+=(Scalar const &rhs) {
  std::visit(
      overload{[&](double &a, complex const &b) { value_ = complex(a) + b; },
               [](auto &&a, auto &&b) { a += b; }},
      value_, rhs.value_);
  return *this;
}

Scalar &Scalar::operator-=(Scalar const &rhs) {
  std::visit(overload{[&](double &a, complex b) { value_ = complex(a) - b; },
                      [](auto &&a, auto &&b) { a -= b; }},
             value_, rhs.value_);
  return *this;
}

Scalar &Scalar::operator*=(Scalar const &rhs) {
  std::visit(overload{[&](double &a, complex b) { value_ = complex(a) * b; },
                      [](auto &&a, auto &&b) { a *= b; }},
             value_, rhs.value_);
  return *this;
}

Scalar &Scalar::operator/=(Scalar const &rhs) {
  std::visit(overload{[&](double &a, complex b) { value_ = complex(a) / b; },
                      [](auto &&a, auto &&b) { a /= b; }},
             value_, rhs.value_);
  return *this;
}

Scalar Scalar::operator-() const {
  return std::visit([](auto &&a) { return Scalar(-a); }, value_);
}

Scalar Scalar::operator+(const Scalar &b) const {
  auto a = *this;
  a += b;
  return a;
}

Scalar Scalar::operator-(const Scalar &b) const {
  auto a = *this;
  a -= b;
  return a;
}

Scalar Scalar::operator*(const Scalar &b) const {
  auto a = *this;
  a *= b;
  return a;
}

Scalar Scalar::operator/(const Scalar &b) const {
  auto a = *this;
  a /= b;
  return a;
}

template <typename T> bool Scalar::is() const {
  return std::holds_alternative<T>(value_);
}
template bool Scalar::is<double>() const;
template bool Scalar::is<complex>() const;

template <> double Scalar::as<double>() const try {
  if (const double *v = std::get_if<double>(&value_)) {
    return *v;
  } else {
    XDIAG_THROW("Cannot convert Scalar holding a value of type "
                "\"complex\" to value of type \"double\"");
  }
} XDIAG_CATCH

template <> complex Scalar::as<complex>() const {
  if (const double *v = std::get_if<double>(&value_)) {
    return complex(*v);
  } else {
    return std::get<complex>(value_);
  }
}
Scalar zero(Scalar s) {
  if (isreal(s)) {
    return Scalar(0.0);
  } else {
    return Scalar(complex(0.0, 0.0));
  }
}

bool isreal(Scalar const &s) { return s.isreal(); }
double real(Scalar const &s) { return s.real(); }
double imag(Scalar const &s) { return s.imag(); }
double abs(Scalar const &s) { return s.abs(); }
Scalar conj(Scalar const &s) { return s.conj(); }

bool isapprox(Scalar const &a, Scalar const &b, double rtol, double atol) {
  return a.isapprox(b, rtol, atol);
}

std::ostream &operator<<(std::ostream &out, Scalar const &v) {
  if (v.isreal()) {
    out << fmt::format("{:.8e}", real(v));
  } else {
    out << fmt::format("({:.8e},{:.8e})", real(v), imag(v));
  }
  return out;
}
std::string to_string(Scalar const &v) { return to_string_generic(v); }

} // namespace xdiag
