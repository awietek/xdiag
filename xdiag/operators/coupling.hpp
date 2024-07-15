#pragma once

#include <string>
#include <variant>
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>

namespace xdiag {

class Coupling {
public:
  using variant_t =
      std::variant<std::string, double, complex, arma::mat, arma::cx_mat>;
  Coupling() = default;
  explicit Coupling(std::string value);
  explicit Coupling(double value);
  explicit Coupling(complex value);
  explicit Coupling(arma::mat const &value);
  explicit Coupling(arma::cx_mat const &value);
  explicit Coupling(variant_t const &value);

  bool isreal() const;
  bool ismatrix() const;
  bool isexplicit() const;

  template <typename coeff_t> bool is() const;
  template <typename coeff_t> coeff_t as() const;

  std::string type() const;

  bool operator==(Coupling const &rhs) const;
  bool operator!=(Coupling const &rhs) const;

  Coupling &operator=(std::string rhs);
  Coupling &operator=(double rhs);
  Coupling &operator=(complex rhs);
  Coupling &operator=(arma::mat const &rhs);
  Coupling &operator=(arma::cx_mat const &rhs);

  // vector space operations
  Coupling &operator+=(Coupling const &rhs);
  Coupling &operator-=(Coupling const &rhs);

  Coupling &operator*=(double rhs);
  Coupling &operator*=(complex rhs);

  variant_t value() const;

private:
  variant_t value_;
};

Coupling &operator/=(Coupling &a, double b);
Coupling &operator/=(Coupling &a, complex b);

Coupling operator+(Coupling const &a, Coupling const &b);
Coupling operator-(Coupling const &a, Coupling const &b);

Coupling operator*(double scalar, Coupling const &x);
Coupling operator*(complex scalar, Coupling const &x);
Coupling operator*(Coupling const &x, double scalar);
Coupling operator*(Coupling const &x, complex scalar);

Coupling operator/(double scalar, Coupling const &x);
Coupling operator/(complex scalar, Coupling const &x);
Coupling operator/(Coupling const &x, double scalar);
Coupling operator/(Coupling const &x, complex scalar);

std::ostream &operator<<(std::ostream &out, Coupling const &cpl);

} // namespace xdiag
