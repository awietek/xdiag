// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <variant>
#include <vector>

#include <xdiag/armadillo.hpp>
#include <xdiag/math/scalar.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Vector is a type-erased column vector: either arma::vec (real) or
// arma::cx_vec (complex). Mixed-type arithmetic widens to complex.
//
// Use is<arma::vec>() / is<arma::cx_vec>() to query; as<T>() to extract
// (as<arma::cx_vec>() always succeeds; as<arma::vec>() throws if complex).
// conj() returns the element-wise complex conjugate as a column vector.
// to_real() narrows back to arma::vec if the imaginary part is zero.
class Vector {
public:
  using value_t = std::variant<arma::vec, arma::cx_vec>;

  Vector() = default;
  Vector(arma::vec const &vec);
  Vector(arma::cx_vec const &vec);
  Vector(std::vector<double> const &vec);
  Vector(std::vector<complex> const &vec);

  bool operator==(Vector const &rhs) const;
  bool operator!=(Vector const &rhs) const;

  // Vector space operations. Mixed real/complex widens the result to complex.
  Vector &operator+=(Vector const &rhs);
  Vector &operator-=(Vector const &rhs);
  Vector &operator*=(Scalar const &rhs);
  Vector &operator/=(Scalar const &rhs);

  Vector operator-() const;
  Vector operator+(Vector const &b) const;
  Vector operator-(Vector const &b) const;
  Vector operator*(Scalar const &rhs) const;
  Vector operator/(Scalar const &rhs) const;

  template <typename T> bool is() const;
  template <typename T> T as() const;

  int64_t size() const;
  bool isreal() const;
  arma::vec real() const;
  arma::vec imag() const; // zero vector for real Vectors
  Vector conj() const;    // element-wise conjugate (column vector)
  Vector to_real(double tol = 1e-12)
      const; // narrows to arma::vec; throws if any |imag| > tol
  bool isapprox(Vector const &y, double rtol = 1e-12,
                double atol = 1e-12) const;

private:
  value_t vec_;
};

bool isreal(Vector const &v);
arma::vec real(Vector const &v);
arma::vec imag(Vector const &v);
Vector conj(Vector const &v);
bool isapprox(Vector const &a, Vector const &b, double rtol = 1e-12,
              double atol = 1e-12);

std::ostream &operator<<(std::ostream &out, Vector const &v);
std::string to_string(Vector const &v);

} // namespace xdiag
