// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "vector.hpp"
#include <xdiag/utils/arma_to_cx.hpp>
#include <xdiag/utils/type_string.hpp>

namespace xdiag {

Vector::Vector(arma::vec const &vec) : vec_(vec) {}
Vector::Vector(arma::cx_vec const &vec) : vec_(vec) {}

template <typename T> bool Vector::is() const {
  return std::holds_alternative<T>(vec_);
}
template bool Vector::is<arma::vec>() const;
template bool Vector::is<arma::cx_vec>() const;

template <> arma::vec Vector::as<arma::vec>() const try {
  if (const arma::vec *m = std::get_if<arma::vec>(&vec_)) {
    return *m;
  } else {
    XDIAG_THROW(fmt::format("Cannot convert Vector holding a value of type "
                            "\"{}\" to value of type \"{}\"",
                            utils::type_string<arma::cx_vec>(),
                            utils::type_string<arma::vec>()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> arma::cx_vec Vector::as<arma::cx_vec>() const {
  if (const arma::vec *m = std::get_if<arma::vec>(&vec_)) {
    return utils::to_cx_vec(*m);
  } else {
    return std::get<arma::cx_vec>(vec_);
  }
}

int64_t Vector::size() const {
  return std::visit([](auto &&v) { return v.size(); }, vec_);
}
bool Vector::isreal() const { return std::holds_alternative<arma::vec>(vec_); }

arma::vec Vector::real() const {
  return std::visit(
      overload{[](arma::vec const &v) { return v; },
               [](arma::cx_vec const &v) { return arma::vec(arma::real(v)); }},
      vec_);
}

arma::vec Vector::imag() const {
  return std::visit(
      overload{[](arma::vec const &v) {
                 return arma::vec(v.n_elem, arma::fill::zeros);
               },
               [](arma::cx_vec const &v) { return arma::vec(arma::imag(v)); }},
      vec_);
}
Vector Vector::hc() const {
  return std::visit(
      overload{[](arma::vec const &v) { return Vector(v); },
               [](arma::cx_vec const &v) { return Vector(arma::conj(v)); }},
      vec_);
}

bool Vector::isapprox(Vector const &y, double rtol, double atol) const {
  return std::visit(overload{
                        [&](arma::vec const &a, arma::vec const &b) {
                          return arma::approx_equal(a, b, "both", atol, rtol);
                        },
                        [&](arma::vec const &a, arma::cx_vec const &&b) {
                          return arma::approx_equal(utils::to_cx_vec(a), b,
                                                    "both", atol, rtol);
                        },
                        [&](arma::cx_vec const &a, arma::vec const &&b) {
                          return arma::approx_equal(a, utils::to_cx_vec(b),
                                                    "both", atol, rtol);
                        },
                        [&](arma::cx_vec const &a, arma::cx_vec const &b) {
                          return arma::approx_equal(a, b, "both", atol, rtol);
                        },
                        [&](auto &&a, auto &&b) { return false; },
                    },
                    vec_, y.vec_);
}

bool Vector::operator==(Vector const &rhs) const {
  return isapprox(rhs, 1e-16, 1e-16);
}
bool Vector::operator!=(Vector const &rhs) const { return !operator==(rhs); }

Vector &Vector::operator+=(Vector const &rhs) {
  std::visit(overload{[&](arma::vec &a, arma::cx_vec b) {
                        vec_ = arma::cx_vec(utils::to_cx_vec(a) + b);
                      },
                      [&](arma::cx_vec &a, arma::vec b) {
                        vec_ = arma::cx_vec(a + utils::to_cx_vec(b));
                      },
                      [](auto &&a, auto &&b) { a += b; }},
             vec_, rhs.vec_);
  return *this;
}
Vector &Vector::operator-=(Vector const &rhs) {
  std::visit(overload{[&](arma::vec &a, arma::cx_vec b) {
                        vec_ = arma::cx_vec(utils::to_cx_vec(a) - b);
                      },
                      [&](arma::cx_vec &a, arma::vec b) {
                        vec_ = arma::cx_vec(a - utils::to_cx_vec(b));
                      },
                      [](auto &&a, auto &&b) { a -= b; }},
             vec_, rhs.vec_);
  return *this;
}

Vector &Vector::operator*=(Scalar const &rhs) {
  if (rhs.isreal()) {
    std::visit(
        overload{[&](arma::vec &m) { vec_ = arma::vec(m * rhs.as<double>()); },
                 [&](arma::cx_vec &m) {
                   vec_ = arma::cx_vec(m * rhs.as<complex>());
                 }},
        vec_);
  } else {
    std::visit(overload{[&](arma::vec &v) {
                          vec_ = arma::cx_vec(v * rhs.as<complex>());
                        },
                        [&](arma::cx_vec &v) {
                          vec_ = arma::cx_vec(v * rhs.as<complex>());
                        }},
               vec_);
  }
  return *this;
}

Vector &Vector::operator/=(Scalar const &rhs) {
  if (rhs.isreal()) {
    std::visit(
        overload{[&](arma::vec &v) { vec_ = arma::vec(v / rhs.as<double>()); },
                 [&](arma::cx_vec &v) {
                   vec_ = arma::cx_vec(v / rhs.as<complex>());
                 }},
        vec_);
  } else {
    std::visit(overload{[&](arma::vec &v) {
                          vec_ = arma::cx_vec(v / rhs.as<complex>());
                        },
                        [&](arma::cx_vec &v) {
                          vec_ = arma::cx_vec(v / rhs.as<complex>());
                        }},
               vec_);
  }
  return *this;
}

Vector Vector::operator-() const {
  return std::visit(
      overload{[](arma::vec const &a) { return Vector(arma::vec(-a)); },
               [](arma::cx_vec const &a) { return Vector(arma::cx_vec(-a)); }},
      vec_);
}

Vector Vector::operator+(Vector const &b) const {
  auto a = *this;
  a += b;
  return a;
}

Vector Vector::operator-(Vector const &b) const {
  auto a = *this;
  a -= b;
  return a;
}

Vector Vector::operator*(Scalar const &b) const {
  auto a = *this;
  a *= b;
  return a;
}

Vector Vector::operator/(Scalar const &b) const {
  auto a = *this;
  a /= b;
  return a;
}

bool isreal(Vector const &s) { return s.isreal(); }
arma::vec real(Vector const &s) { return s.real(); }
arma::vec imag(Vector const &s) { return s.imag(); }
Vector hc(Vector const &s) { return s.hc(); }

bool isapprox(Vector const &a, Vector const &b, double rtol, double atol) {
  return a.isapprox(b, rtol, atol);
}

std::ostream &operator<<(std::ostream &out, Vector const &vec) {
  if (vec.isreal()) {
    vec.as<arma::vec>().brief_print(out);
  } else {
    vec.as<arma::cx_vec>().brief_print(out);
  }
  return out;
}
std::string to_string(Vector const &vec) { return to_string_generic(vec); }

} // namespace xdiag
