#include "matrix.hpp"
#include <xdiag/utils/type_string.hpp>
#include <xdiag/utils/arma_to_cx.hpp>

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
                            utils::type_string<arma::cx_mat>(),
                            utils::type_string<arma::mat>()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> arma::cx_mat Matrix::as<arma::cx_mat>() const {
  if (const arma::mat *m = std::get_if<arma::mat>(&mat_)) {
    return to_cx_mat(*m);
  } else {
    return std::get<arma::cx_mat>(mat_);
  }
}

bool Matrix::isreal() const { return std::holds_alternative<arma::mat>(mat_); }

arma::mat Matrix::real() const {
  return std::visit(
      overload{[](arma::mat const &m) { return m; },
               [](arma::cx_mat const &m) { return arma::mat(arma::real(m)); }},
      mat_);
}

arma::mat Matrix::imag() const {
  return std::visit(
      overload{[](arma::mat const &m) {
                 return arma::mat(m.n_rows, m.n_cols, arma::fill::zeros);
               },
               [](arma::cx_mat const &m) { return arma::mat(arma::imag(m)); }},
      mat_);
}
Matrix Matrix::hc() const {
  return std::visit(overload{[](arma::mat const &m) {
                               return Matrix(arma::mat(arma::trans(m)));
                             },
                             [](arma::cx_mat const &m) {
                               return Matrix(arma::cx_mat(arma::trans(m)));
                             }},
                    mat_);
}

bool Matrix::isapprox(Matrix const &y, double rtol, double atol) const {
  return std::visit(
      overload{
          [&](arma::mat const &a, arma::mat const &b) {
            return arma::approx_equal(a, b, "both", atol, rtol);
          },
          [&](arma::mat const &a, arma::cx_mat const &&b) {
            return arma::approx_equal(to_cx_mat(a), b, "both", atol, rtol);
          },
          [&](arma::cx_mat const &a, arma::mat const &&b) {
            return arma::approx_equal(a, to_cx_mat(b), "both", atol, rtol);
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

Matrix &Matrix::operator+=(Matrix const &rhs) {
  std::visit(overload{[&](arma::mat &a, arma::cx_mat b) {
                        mat_ = arma::cx_mat(to_cx_mat(a) + b);
                      },
                      [&](arma::cx_mat &a, arma::mat b) {
                        mat_ = arma::cx_mat(a + to_cx_mat(b));
                      },
                      [](auto &&a, auto &&b) { a += b; }},
             mat_, rhs.mat_);
  return *this;
}
Matrix &Matrix::operator-=(Matrix const &rhs) {
  std::visit(overload{[&](arma::mat &a, arma::cx_mat b) {
                        mat_ = arma::cx_mat(to_cx_mat(a) - b);
                      },
                      [&](arma::cx_mat &a, arma::mat b) {
                        mat_ = arma::cx_mat(a - to_cx_mat(b));
                      },
                      [](auto &&a, auto &&b) { a -= b; }},
             mat_, rhs.mat_);
  return *this;
}

Matrix &Matrix::operator*=(Scalar const &rhs) {
  if (rhs.isreal()) {
    std::visit(
        overload{[&](arma::mat &m) { mat_ = arma::mat(m * rhs.as<double>()); },
                 [&](arma::cx_mat &m) {
                   mat_ = arma::cx_mat(m * rhs.as<complex>());
                 }},
        mat_);
  } else {
    std::visit(overload{[&](arma::mat &m) {
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
    std::visit(
        overload{[&](arma::mat &m) { mat_ = arma::mat(m / rhs.as<double>()); },
                 [&](arma::cx_mat &m) {
                   mat_ = arma::cx_mat(m / rhs.as<complex>());
                 }},
        mat_);
  } else {
    std::visit(overload{[&](arma::mat &m) {
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
      overload{[](arma::mat const &a) { return Matrix(arma::mat(-a)); },
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

bool isreal(Matrix const &s) { return s.isreal(); }
arma::mat real(Matrix const &s) { return s.real(); }
arma::mat imag(Matrix const &s) { return s.imag(); }
Matrix hc(Matrix const &s) { return s.hc(); }

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
