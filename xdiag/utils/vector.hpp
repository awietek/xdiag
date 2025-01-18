#pragma once

#include <string>
#include <variant>
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/utils/scalar.hpp>

namespace xdiag {

class Vector {
public:
  using value_t = std::variant<arma::vec, arma::cx_vec>;

  Vector() = default;
  Vector(arma::vec const &vec);
  Vector(arma::cx_vec const &vec);

  bool operator==(Vector const &rhs) const;
  bool operator!=(Vector const &rhs) const;

  // vector space operations
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

  bool isreal() const;
  arma::vec real() const;
  arma::vec imag() const;
  Vector hc() const;
  bool isapprox(Vector const &y, double rtol = 1e-12,
                double atol = 1e-12) const;

private:
  value_t vec_;
};

bool isreal(Vector const &s);
arma::vec real(Vector const &s);
arma::vec imag(Vector const &s);
Vector hc(Vector const &s);
bool isapprox(Vector const &a, Vector const &b, double rtol = 1e-10,
              double atol = 0);

std::ostream &operator<<(std::ostream &out, Vector const &cpl);
std::string to_string(Vector const &cpl);

} // namespace xdiag
