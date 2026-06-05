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

  XDIAG_API Vector() = default;
  XDIAG_API Vector(arma::vec const &vec);
  XDIAG_API Vector(arma::cx_vec const &vec);
  XDIAG_API Vector(std::vector<double> const &vec);
  XDIAG_API Vector(std::vector<complex> const &vec);

  XDIAG_API bool operator==(Vector const &rhs) const;
  XDIAG_API bool operator!=(Vector const &rhs) const;

  // Vector space operations. Mixed real/complex widens the result to complex.
  XDIAG_API Vector &operator+=(Vector const &rhs);
  XDIAG_API Vector &operator-=(Vector const &rhs);
  XDIAG_API Vector &operator*=(Scalar const &rhs);
  XDIAG_API Vector &operator/=(Scalar const &rhs);

  XDIAG_API Vector operator-() const;
  XDIAG_API Vector operator+(Vector const &b) const;
  XDIAG_API Vector operator-(Vector const &b) const;
  XDIAG_API Vector operator*(Scalar const &rhs) const;
  XDIAG_API Vector operator/(Scalar const &rhs) const;

  template <typename T> XDIAG_API bool is() const;
  template <typename T> XDIAG_API T as() const;

  XDIAG_API int64_t size() const;
  XDIAG_API bool isreal() const;
  XDIAG_API arma::vec real() const;
  XDIAG_API arma::vec imag() const; // zero vector for real Vectors
  XDIAG_API Vector conj() const;    // element-wise conjugate (column vector)
  XDIAG_API Vector to_real(double tol = 1e-12)
      const; // narrows to arma::vec; throws if any |imag| > tol
  XDIAG_API bool isapprox(Vector const &y, double rtol = 1e-12,
                          double atol = 1e-12) const;

private:
  value_t vec_;
};

XDIAG_API bool isreal(Vector const &v);
XDIAG_API arma::vec real(Vector const &v);
XDIAG_API arma::vec imag(Vector const &v);
XDIAG_API Vector conj(Vector const &v);
XDIAG_API bool isapprox(Vector const &a, Vector const &b, double rtol = 1e-12,
                        double atol = 1e-12);

XDIAG_API std::ostream &operator<<(std::ostream &out, Vector const &v);
XDIAG_API std::string to_string(Vector const &v);

} // namespace xdiag
