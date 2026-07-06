// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <variant>

#include <xdiag/armadillo.hpp>
#include <xdiag/math/scalar.hpp>
#include <xdiag/math/vector.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Matrix is a type-erased dense matrix: either arma::mat (real) or
// arma::cx_mat (complex). Mixed-type arithmetic widens to complex.
//
// Use is<arma::mat>() / is<arma::cx_mat>() to query; as<T>() to extract
// (as<arma::cx_mat>() always succeeds; as<arma::mat>() throws if complex).
// hc() returns the conjugate transpose (Hermitian conjugate).
// to_real() narrows back to arma::mat if the imaginary part is zero.
class Matrix {
public:
  using value_t = std::variant<arma::mat, arma::cx_mat>;

  Matrix() = default;
  Matrix(arma::mat const &mat);
  Matrix(arma::cx_mat const &mat);

  bool operator==(Matrix const &rhs) const;
  bool operator!=(Matrix const &rhs) const;
  bool operator<(Matrix const &rhs) const;

  // Linear combination (scalar multiply/add). Mixed real/complex widens to
  // complex.
  Matrix &operator+=(Matrix const &rhs);
  Matrix &operator-=(Matrix const &rhs);
  Matrix &operator*=(Scalar const &rhs);
  Matrix &operator/=(Scalar const &rhs);

  Matrix operator-() const;
  Matrix operator+(Matrix const &b) const;
  Matrix operator-(Matrix const &b) const;
  Matrix operator*(Scalar const &rhs) const;
  Matrix operator/(Scalar const &rhs) const;

  // Matrix-matrix and matrix-vector products. Widens to complex if either side
  // is complex.
  Matrix operator*(Matrix const &rhs) const;
  Vector operator*(Vector const &rhs) const;

  template <typename T> bool is() const;
  template <typename T> T as() const;

  int64_t n_rows() const;
  int64_t n_cols() const;
  bool isreal() const;
  arma::mat real() const;
  arma::mat imag() const; // zero matrix for real Matrices
  Matrix hc() const;      // conjugate transpose
  Matrix to_real(double tol = 1e-12)
      const; // narrows to arma::mat; throws if any |imag| > tol
  bool isapprox(Matrix const &y, double rtol = 1e-12,
                double atol = 1e-12) const;

private:
  value_t mat_;
};

bool isreal(Matrix const &m);
arma::mat real(Matrix const &m);
arma::mat imag(Matrix const &m);
Matrix hc(Matrix const &m);
bool isapprox(Matrix const &a, Matrix const &b, double rtol = 1e-12,
              double atol = 1e-12);

std::ostream &operator<<(std::ostream &out, Matrix const &mat);
std::string to_string(Matrix const &mat);

} // namespace xdiag
