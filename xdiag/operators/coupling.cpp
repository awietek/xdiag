#include "coupling.hpp"

namespace xdiag {

Coupling::Coupling(std::string value) : value_(value) {}
Coupling::Coupling(double value) : value_(value) {}
Coupling::Coupling(complex value) : value_(value) {}
Coupling::Coupling(arma::mat const &value) : value_(value) {}
Coupling::Coupling(arma::cx_mat const &value) : value_(value) {}
Coupling::Coupling(variant_t const &value) : value_(value) {}

bool Coupling::isreal() const try {
  return std::visit(overload{
                        [](std::string const &) {
                          XDIAG_THROW("Cannot determine whether coupling is "
                                      "real, since it is of std::string type");
                          return false;
                        },
                        [](double const &) { return true; },
                        [](complex const &) { return false; },
                        [](arma::mat const &) { return true; },
                        [](arma::cx_mat const &) { return false; },
                    },
                    value_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool Coupling::ismatrix() const {
  return std::visit(overload{
                        [](std::string const &) { return false; },
                        [](double const &) { return false; },
                        [](complex const &) { return false; },
                        [](arma::mat const &) { return true; },
                        [](arma::cx_mat const &) { return true; },
                    },
                    value_);
}

bool Coupling::isexplicit() const {
  return std::visit(overload{
                        [](std::string const &) { return false; },
                        [](double const &) { return true; },
                        [](complex const &) { return true; },
                        [](arma::mat const &) { return true; },
                        [](arma::cx_mat const &) { return true; },
                    },
                    value_);
}

template <typename coeff_t> bool Coupling::is() const {
  return std::holds_alternative<coeff_t>(value_);
}
template bool Coupling::is<std::string>() const;
template bool Coupling::is<double>() const;
template bool Coupling::is<complex>() const;
template bool Coupling::is<arma::mat>() const;
template bool Coupling::is<arma::cx_mat>() const;

template <typename coeff_t> coeff_t Coupling::as() const try {
  if (const coeff_t *pval = std::get_if<coeff_t>(&value_)) {
    return *pval;
  } else {
    XDIAG_THROW(fmt::format(
        "Cannot retrieve desired value, since Coupling is of type \"{}\"",
        type()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template std::string Coupling::as<std::string>() const;
template double Coupling::as<double>() const;
template arma::mat Coupling::as<arma::mat>() const;

template <> complex Coupling::as<complex>() const try {
  if (const double *pval = std::get_if<double>(&value_)) {
    return (complex)*pval;
  } else if (const complex *pval = std::get_if<complex>(&value_)) {
    return *pval;
  } else {
    XDIAG_THROW(fmt::format(
        "Cannot retrieve desired value, since Coupling is of type \"{}\"",
        type()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> arma::cx_mat Coupling::as<arma::cx_mat>() const try {
  if (const arma::mat *pval = std::get_if<arma::mat>(&value_)) {
    return to_cx_mat(*pval);
  } else if (const arma::cx_mat *pval = std::get_if<arma::cx_mat>(&value_)) {
    return *pval;
  } else {
    XDIAG_THROW(fmt::format(
        "Cannot retrieve desired value, since Coupling is of type \"{}\"",
        type()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::string Coupling::type() const {
  return std::visit(overload{
                        [](std::string const &) { return "std::string"; },
                        [](double const &) { return "double"; },
                        [](complex const &) { return "complex"; },
                        [](arma::mat const &) { return "arma::mat"; },
                        [](arma::cx_mat const &) { return "arma::cx_mat"; },
                    },
                    value_);
}

bool Coupling::operator==(Coupling const &rhs) const {
  return std::visit(
      overload{
          [](std::string const &a, std::string const &b) { return a == b; },
          [](double const &a, double const &b) { return a == b; },
          [](complex const &a, complex const &b) { return a == b; },
          [](arma::mat const &a, arma::mat const &b) {
            return arma::all(vectorise(a) == vectorise(b));
          },
          [](arma::cx_mat const &a, arma::cx_mat const &b) {
            return arma::all(vectorise(a) == vectorise(b));
          },
          [](auto const &, auto const &) { return false; }},
      value_, rhs.value_);
}

bool Coupling::operator!=(Coupling const &rhs) const {
  return !operator==(rhs);
}

Coupling &Coupling::operator=(std::string rhs) {
  value_ = rhs;
  return *this;
}
Coupling &Coupling::operator=(double rhs) {
  value_ = rhs;
  return *this;
}
Coupling &Coupling::operator=(complex rhs) {
  value_ = rhs;
  return *this;
}
Coupling &Coupling::operator=(arma::mat const &rhs) {
  value_ = rhs;
  return *this;
}
Coupling &Coupling::operator=(arma::cx_mat const &rhs) {
  value_ = rhs;
  return *this;
}

Coupling &Coupling::operator+=(Coupling const &rhs) try {
  std::visit(
      overload{[](double &a, double const &b) { a += b; },
               [&](double &a, complex const &b) { value_ = (complex)a + b; },
               [](complex &a, double const &b) { a += b; },
               [](complex &a, complex const &b) { a += b; },
               [](arma::mat &a, arma::mat const &b) { a += b; },
               [&](arma::mat &a, arma::cx_mat const &b) {
                 auto acplx = to_cx_mat(a);
                 value_ = arma::cx_mat(acplx + b);
               },
               [](arma::cx_mat &a, arma::mat const &b) {
                 auto bcplx = to_cx_mat(b);
                 a += bcplx;
               },
               [](arma::cx_mat &a, arma::cx_mat const &b) { a += b; },
               [&](auto &&, auto &&) {
                 XDIAG_THROW(fmt::format(
                     "Cannot add couplings of type \"{}\" and \"{}\"", type(),
                     rhs.type()));
               }},
      value_, rhs.value_);
  return *this;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling &Coupling::operator-=(Coupling const &rhs) try {
  std::visit(
      overload{[](double &a, double const &b) { a -= b; },
               [&](double &a, complex const &b) { value_ = (complex)a - b; },
               [](complex &a, double const &b) { a -= b; },
               [](complex &a, complex const &b) { a -= b; },
               [](arma::mat &a, arma::mat const &b) { a -= b; },
               [&](arma::mat &a, arma::cx_mat const &b) {
                 auto acplx = to_cx_mat(a);
                 value_ = arma::cx_mat(acplx - b);
               },
               [](arma::cx_mat &a, arma::mat const &b) {
                 auto bcplx = to_cx_mat(b);
                 a -= bcplx;
               },
               [](arma::cx_mat &a, arma::cx_mat const &b) { a -= b; },
               [&](auto &&, auto &&) {
                 XDIAG_THROW(fmt::format(
                     "Cannot subtract couplings of type \"{}\" and \"{}\"",
                     type(), rhs.type()));
               }},
      value_, rhs.value_);
  return *this;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling &Coupling::operator*=(double rhs) try {
  std::visit(
      overload{
          [](std::string const &) {
            XDIAG_THROW(
                "Cannot scalar multiply a coupling of type \"std::string\"");
          },
          [&](double &a) { a *= rhs; },
          [&](complex &a) { a *= rhs; },
          [&](arma::mat &a) { a *= rhs; },
          [&](arma::cx_mat &a) { a *= rhs; },
      },
      value_);
  return *this;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling &Coupling::operator*=(complex rhs) try {
  std::visit(
      overload{
          [](std::string const &) {
            XDIAG_THROW(
                "Cannot scalar multiply a coupling of type \"std::string\"");
          },
          [&](double &a) { value_ = (complex)a * rhs; },
          [&](complex &a) { a *= rhs; },
          [&](arma::mat &a) { value_ = arma::cx_mat(to_cx_mat(a) * rhs); },
          [&](arma::cx_mat &a) { a *= rhs; },
      },
      value_);
  return *this;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling::variant_t Coupling::value() const { return value_; }

Coupling &operator/=(Coupling &a, double b) try {
  return a *= (1 / b);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling &operator/=(Coupling &a, complex b) try {
  return a *= (1 / b);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling operator+(Coupling const &a, Coupling const &b) try {
  Coupling ret = a;
  ret += b;
  return ret;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling operator-(Coupling const &a, Coupling const &b) try {
  Coupling ret = a;
  ret -= b;
  return ret;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling operator*(double scalar, Coupling const &x) try {
  Coupling ret = x;
  ret *= scalar;
  return ret;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling operator*(complex scalar, Coupling const &x) try {
  Coupling ret = x;
  ret *= scalar;
  return ret;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling operator*(Coupling const &x, double scalar) try {
  Coupling ret = x;
  ret *= scalar;
  return ret;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling operator*(Coupling const &x, complex scalar) try {
  Coupling ret = x;
  ret *= scalar;
  return ret;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling operator/(double scalar, Coupling const &x) try {
  Coupling ret = x;
  ret /= scalar;
  return ret;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling operator/(complex scalar, Coupling const &x) try {
  Coupling ret = x;
  ret /= scalar;
  return ret;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling operator/(Coupling const &x, double scalar) try {
  Coupling ret = x;
  ret /= scalar;
  return ret;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Coupling operator/(Coupling const &x, complex scalar) try {
  Coupling ret = x;
  ret /= scalar;
  return ret;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::ostream &operator<<(std::ostream &out, Coupling const &cpl) {
  std::visit([&out](auto &&val) { out << val; }, cpl.value());
  return out;
}
std::string to_string(Coupling const &cpl) { return to_string_generic(cpl); }
} // namespace xdiag
