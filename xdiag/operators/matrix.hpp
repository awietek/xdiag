#pragma once

#include <string>
#include <variant>
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/scalar.hpp>

namespace xdiag {

class Matrix {
public:
  using value_t = std::variant<arma::mat, arma::cx_mat>;

  Matrix() = default;
  explicit Matrix(arma::mat const &mat);
  explicit Matrix(arma::cx_mat const &mat);

  bool operator==(Matrix const &rhs) const;
  bool operator!=(Matrix const &rhs) const;

  // vector space operations
  Matrix &operator+=(Matrix const &rhs);
  Matrix &operator-=(Matrix const &rhs);
  Matrix &operator*=(Scalar const &rhs);
  Matrix &operator/=(Scalar const &rhs);

  Matrix operator-() const;
  Matrix operator+(Matrix const &b) const;
  Matrix operator-(Matrix const &b) const;
  Matrix &operator*(Scalar const &rhs);
  Matrix &operator/(Scalar const &rhs);

  template <typename T> bool is() const;
  template <typename T> T as() const;

  bool isreal() const;
  arma::mat real() const;
  arma::mat imag() const;
  Matrix hc() const;
  bool isapprox(Matrix const &y, double rtol = 1e-12,
                double atol = 1e-12) const;

private:
  value_t mat_;
};

bool isreal(Matrix const &s);
arma::mat real(Matrix const &s);
arma::mat imag(Matrix const &s);
Matrix hc(Matrix const &s);
bool isapprox(Matrix const &a, Matrix const &b, double rtol = 1e-10,
              double atol = 0);

std::ostream &operator<<(std::ostream &out, Matrix const &cpl);
std::string to_string(Matrix const &cpl);

} // namespace xdiag
