#include "matrix.hpp"
#include <xdiag/utils/type_string.hpp>

namespace xdiag {

Matrix::Matrix(arma::mat const &mat) : mat_(mat) {}
Matrix::Matrix(arma::cx_mat const &mat) : mat_(mat) {}

template <typename T> bool Matrix::is() const {
  return std::holds_alternative<T>(mat_);
}
template bool Matrix::is<arma::mat>() const;
template bool Matrix::is<arma::cx_mat>() const;

template arma::mat Matrix::as<arma::mat>() const try {
  if (const arma::mat *m = std::get_if<arma::mat>(&mat_)) {
    return *m;
  } else {
    XDIAG_THROW(fmt::format("Cannot convert Matrix holding a value of type "
                            "\"{}\" to value of type \"{}\"",
                            type_string<arma::cx_mat>(),
                            type_string<arma::mat>()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

static arma::cx_mat to_cx_mat(arma::mat const &mat) {
  return arma::cx_mat(mat,
                      arma::mat(mat.n_rows, mat.n_cols, arma::fill::zeros));
}

template arma::cx_mat Matrix::as<arma::cx_mat>() const {
  if (const arma::mat *m = std::get_if<arma::mat>(&mat_)) {
    return to_cx_mat(*m);
  } else {
    return std::get<arma::cx_mat>(mat_)
  }
}

bool Matrix::isreal() const { return std::holds_alternative<arma::mat>(mat_); }

arma::mat Matrix::real() const {
  return std::visit(
      overload{[](arma::mat const &m) { return m; },
               [](arma::cx_mat const &m) { return arma::real(m); }},
      mat_);
}

arma::mat Matrix::imag() const {
  return std::visit(
      overload{[](arma::mat const &m) {
                 return arma::mat(m.n_rows, m.n_cols, arma::fill::zeros);
               },
               [](arma::cx_mat const &m) { return arma::imag(m); }},
      mat_);
}
Matrix Matrix::hc() const {
  return Matrix(std::visit([](auto &&m) { return arma::trans(m); }, mat_));
}

bool Matrix::isapprox(Matrix const &y, double rtol = 1e-12,
                      double atol = 1e-12) const {
  return std::visit([](auto &&a, auto &&b) {
    return arma::approx_equal(a, b, "both", atol, rtol);
  });
}

bool Matrix::operator==(Matrix const &rhs) const {
  return isapprox(mat_, rhs.mat_, 1e-16, 1e-16);
}
bool Matrix::operator!=(Matrix const &rhs) const { return !operator==(rhs); }

Matrix &Matrix::operator+=(Matrix const &rhs) {
  std::visit(
      overload{[&](arma::mat &a, arma::cx_mat b) { mat_ = to_cx_mat(a) + b; },
               [&](arma::cx_mat &a, arma::mat b) { mat_ = a + to_cx_mat(b); },
               [](auto &&a, auto &&b) { a += b; }},
      mat_, rhs.mat_);
  return *this;
}
Matrix &Matrix::operator-=(Matrix const &rhs) {
  std::visit(
      overload{[&](arma::mat &a, arma::cx_mat b) { mat_ = to_cx_mat(a) - b; },
               [&](arma::cx_mat &a, arma::mat b) { mat_ = a - to_cx_mat(b); },
               [](auto &&a, auto &&b) { a -= b; }},
      mat_, rhs.mat_);
  return *this;
}

Matrix &Matrix::operator*=(Scalar const &rhs) {
  if (rhs.isreal()) {
    std::visit([&](auto &&m) { m *= rhs.as<double>(); });
  } else {
    std::visit(
        overload{[&](arma::mat &m) { mat_ = to_cx_mat(m) * rhs.as<complex>(); },
                 [&](arma::cx_mat &m) { m *= rhs.as<complex>(); }},
        mat_);
  }
  return *this;
}
Matrix &Matrix::operator/=(Scalar const &rhs) {
  if (rhs.isreal()) {
    std::visit([&](auto &&m) { m /= rhs.as<double>(); });
  } else {
    std::visit(
        overload{[&](arma::mat &m) { mat_ = to_cx_mat(m) / rhs.as<complex>(); },
                 [&](arma::cx_mat &m) { m /= rhs.as<complex>(); }},
        mat_);
  }
  return *this;
}

Matrix Matrix::operator-() const {
  return std::visit([](auto &&a) { return -a; }, mat_);
}

Matrix Matrix::operator+(const Matrix &b) const {
  auto a = *this;
  a += b;
  return a;
}

Matrix Matrix::operator-(const Matrix &b) const {
  auto a = *this;
  a -= b;
  return a;
}

Matrix Matrix::operator*(const Scalar &b) const {
  auto a = *this;
  a *= b;
  return a;
}

Matrix Matrix::operator/(const Scalar &b) const {
  auto a = *this;
  a /= b;
  return a;
}

bool isreal(Matrix const &s) { return s.isreal(); }
bool real(Matrix const &s) { return s.real(); }
bool imag(Matrix const &s) { return s.imag(); }
Matrix hc(Matrix const &s) { return s.hc(); }

bool isapprox(Matrix const &a, Matrix const &b, double rtol, double atol) {
  a.isapprox(b, rtol, atol);
}

std::ostream &operator<<(std::ostream &out, Matrix const &mat) {
  if (mat.isreal()) {
    mat.as<arma::mat>.brief_print(out);
  } else {
    mat.as<arma::cx_mat>.brief_print(out);
  }
}
std::string to_string(Matrix const &mat) { return to_string_generic(mat); }

} // namespace xdiag
