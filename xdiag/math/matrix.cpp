// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "matrix.hpp"
#include <algorithm>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>
#include <xdiag/utils/type_name.hpp>
#include <xdiag/utils/variants.hpp>

namespace xdiag {

Matrix::Matrix(arma::mat const &mat) : mat_(mat) {}
Matrix::Matrix(arma::cx_mat const &mat) : mat_(mat) {}

template <typename T> bool Matrix::is() const {
  return std::holds_alternative<T>(mat_);
}
template bool Matrix::is<arma::mat>() const;
template bool Matrix::is<arma::cx_mat>() const;

template <> arma::mat Matrix::as<arma::mat>() const try {
  if (const arma::mat *m = std::get_if<arma::mat>(&mat_)) {
    return *m;
  } else {
    XDIAG_THROW(fmt::format("Cannot convert Matrix holding a value of type "
                            "\"{}\" to value of type \"{}\"",
                            utils::get_type_name<arma::cx_mat>(),
                            utils::get_type_name<arma::mat>()));
  }
}
XDIAG_CATCH

template <> arma::cx_mat Matrix::as<arma::cx_mat>() const {
  if (const arma::mat *m = std::get_if<arma::mat>(&mat_)) {
    return arma::conv_to<arma::cx_mat>::from(*m);
  } else {
    return std::get<arma::cx_mat>(mat_);
  }
}

bool Matrix::isreal() const { return std::holds_alternative<arma::mat>(mat_); }

int64_t Matrix::n_rows() const {
  return std::visit([](auto &&m) { return (int64_t)m.n_rows; }, mat_);
}
int64_t Matrix::n_cols() const {
  return std::visit([](auto &&m) { return (int64_t)m.n_cols; }, mat_);
}

arma::mat Matrix::real() const {
  return std::visit(utils::overload{[](arma::mat const &m) { return m; },
                                    [](arma::cx_mat const &m) {
                                      return arma::mat(arma::real(m));
                                    }},
                    mat_);
}

arma::mat Matrix::imag() const {
  return std::visit(utils::overload{[](arma::mat const &m) {
                                      return arma::mat(m.n_rows, m.n_cols,
                                                       arma::fill::zeros);
                                    },
                                    [](arma::cx_mat const &m) {
                                      return arma::mat(arma::imag(m));
                                    }},
                    mat_);
}
Matrix Matrix::hc() const {
  return std::visit(
      utils::overload{
          [](arma::mat const &m) { return Matrix(arma::mat(arma::trans(m))); },
          [](arma::cx_mat const &m) {
            return Matrix(arma::cx_mat(arma::trans(m)));
          }},
      mat_);
}

Matrix Matrix::operator*(Matrix const &rhs) const {
  if (isreal() && rhs.isreal()) {
    return Matrix(arma::mat(as<arma::mat>() * rhs.as<arma::mat>()));
  } else {
    return Matrix(arma::cx_mat(as<arma::cx_mat>() * rhs.as<arma::cx_mat>()));
  }
}

Vector Matrix::operator*(Vector const &rhs) const {
  if (isreal() && rhs.isreal()) {
    return Vector(arma::vec(as<arma::mat>() * rhs.as<arma::vec>()));
  } else {
    return Vector(arma::cx_vec(as<arma::cx_mat>() * rhs.as<arma::cx_vec>()));
  }
}

Matrix Matrix::to_real(double tol) const try {
  if (isreal()) {
    return *this;
  }
  arma::cx_mat const &m = std::get<arma::cx_mat>(mat_);
  if (arma::norm(arma::imag(m), "inf") > tol) {
    XDIAG_THROW(
        fmt::format("Cannot narrow Matrix to real: max imaginary element {} "
                    "exceeds tolerance {}",
                    arma::norm(arma::imag(m), "inf"), tol));
  }
  return Matrix(arma::mat(arma::real(m)));
}
XDIAG_CATCH

bool Matrix::isapprox(Matrix const &y, double rtol, double atol) const {
  return std::visit(
      utils::overload{
          [&](arma::mat const &a, arma::mat const &b) {
            return arma::approx_equal(a, b, "both", atol, rtol);
          },
          [&](arma::mat const &a, arma::cx_mat const &b) {
            return arma::approx_equal(arma::conv_to<arma::cx_mat>::from(a), b,
                                      "both", atol, rtol);
          },
          [&](arma::cx_mat const &a, arma::mat const &b) {
            return arma::approx_equal(a, arma::conv_to<arma::cx_mat>::from(b),
                                      "both", atol, rtol);
          },
          [&](arma::cx_mat const &a, arma::cx_mat const &b) {
            return arma::approx_equal(a, b, "both", atol, rtol);
          },
          [&](auto &&a, auto &&b) { return false; },
      },
      mat_, y.mat_);
}

bool Matrix::operator==(Matrix const &rhs) const {
  return isapprox(rhs, 1e-16, 1e-16);
}
bool Matrix::operator!=(Matrix const &rhs) const { return !operator==(rhs); }

// Lexicographic ordering: real < complex; within same type, compare
// element-by-element in column-major (Armadillo memory) order.
bool Matrix::operator<(Matrix const &rhs) const {
  // real matrices come before complex
  if (isreal() != rhs.isreal())
    return isreal() >
           rhs.isreal(); // real (true) > complex (false) → real first
  // Compare by shape first
  if (n_rows() != rhs.n_rows())
    return n_rows() < rhs.n_rows();
  if (n_cols() != rhs.n_cols())
    return n_cols() < rhs.n_cols();
  // Same shape: lexicographic comparison of elements (column-major memory)
  if (isreal()) {
    auto const &a = as<arma::mat>();
    auto const &b = rhs.as<arma::mat>();
    int64_t s = n_rows() * n_cols();
    return std::lexicographical_compare(a.memptr(), a.memptr() + s, b.memptr(),
                                        b.memptr() + s);
  } else {
    auto const &a = as<arma::cx_mat>();
    auto const &b = rhs.as<arma::cx_mat>();
    int64_t s = n_rows() * n_cols() * 2;
    auto pa = reinterpret_cast<const double *>(a.memptr());
    auto pb = reinterpret_cast<const double *>(b.memptr());
    return std::lexicographical_compare(pa, pa + s, pb, pb + s);
  }
}

Matrix &Matrix::operator+=(Matrix const &rhs) {
  std::visit(utils::overload{[&](arma::mat &a, arma::cx_mat b) {
                               mat_ = arma::cx_mat(
                                   arma::conv_to<arma::cx_mat>::from(a) + b);
                             },
                             [&](arma::cx_mat &a, arma::mat b) {
                               mat_ = arma::cx_mat(
                                   a + arma::conv_to<arma::cx_mat>::from(b));
                             },
                             [](auto &&a, auto &&b) { a += b; }},
             mat_, rhs.mat_);
  return *this;
}
Matrix &Matrix::operator-=(Matrix const &rhs) {
  std::visit(utils::overload{[&](arma::mat &a, arma::cx_mat b) {
                               mat_ = arma::cx_mat(
                                   arma::conv_to<arma::cx_mat>::from(a) - b);
                             },
                             [&](arma::cx_mat &a, arma::mat b) {
                               mat_ = arma::cx_mat(
                                   a - arma::conv_to<arma::cx_mat>::from(b));
                             },
                             [](auto &&a, auto &&b) { a -= b; }},
             mat_, rhs.mat_);
  return *this;
}

Matrix &Matrix::operator*=(Scalar const &rhs) {
  if (rhs.isreal()) {
    std::visit(utils::overload{[&](arma::mat &m) {
                                 mat_ = arma::mat(m * rhs.as<double>());
                               },
                               [&](arma::cx_mat &m) {
                                 mat_ = arma::cx_mat(m * rhs.as<complex>());
                               }},
               mat_);
  } else {
    std::visit(utils::overload{[&](arma::mat &m) {
                                 mat_ = arma::cx_mat(m * rhs.as<complex>());
                               },
                               [&](arma::cx_mat &m) {
                                 mat_ = arma::cx_mat(m * rhs.as<complex>());
                               }},
               mat_);
  }
  return *this;
}
Matrix &Matrix::operator/=(Scalar const &rhs) {
  if (rhs.isreal()) {
    std::visit(utils::overload{[&](arma::mat &m) {
                                 mat_ = arma::mat(m / rhs.as<double>());
                               },
                               [&](arma::cx_mat &m) {
                                 mat_ = arma::cx_mat(m / rhs.as<complex>());
                               }},
               mat_);
  } else {
    std::visit(utils::overload{[&](arma::mat &m) {
                                 mat_ = arma::cx_mat(m / rhs.as<complex>());
                               },
                               [&](arma::cx_mat &m) {
                                 mat_ = arma::cx_mat(m / rhs.as<complex>());
                               }},
               mat_);
  }
  return *this;
}

Matrix Matrix::operator-() const {
  return std::visit(
      utils::overload{
          [](arma::mat const &a) { return Matrix(arma::mat(-a)); },
          [](arma::cx_mat const &a) { return Matrix(arma::cx_mat(-a)); }},
      mat_);
}

Matrix Matrix::operator+(Matrix const &b) const {
  auto a = *this;
  a += b;
  return a;
}

Matrix Matrix::operator-(Matrix const &b) const {
  auto a = *this;
  a -= b;
  return a;
}

Matrix Matrix::operator*(Scalar const &b) const {
  auto a = *this;
  a *= b;
  return a;
}

Matrix Matrix::operator/(Scalar const &b) const {
  auto a = *this;
  a /= b;
  return a;
}

bool isreal(Matrix const &m) { return m.isreal(); }
arma::mat real(Matrix const &m) { return m.real(); }
arma::mat imag(Matrix const &m) { return m.imag(); }
Matrix hc(Matrix const &m) { return m.hc(); }

bool isapprox(Matrix const &a, Matrix const &b, double rtol, double atol) {
  return a.isapprox(b, rtol, atol);
}

std::ostream &operator<<(std::ostream &out, Matrix const &mat) {
  if (mat.isreal()) {
    mat.as<arma::mat>().brief_print(out);
  } else {
    mat.as<arma::cx_mat>().brief_print(out);
  }
  return out;
}
std::string to_string(Matrix const &mat) { return to_string_generic(mat); }

} // namespace xdiag
