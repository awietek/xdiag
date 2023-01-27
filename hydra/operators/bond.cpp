#include "bond.h"

#include <algorithm>
#include <hydra/common.h>
#include <hydra/utils/close.h>
#include <set>

namespace hydra {

// Constructors with type name
Bond::Bond(std::string type, int site)
    : type_(type), matrix_(), coupling_(1.0),
      coupling_name_("HYDRA_COUPLING_NAMELESS"), sites_({site}) {}

Bond::Bond(std::string type, std::vector<int> const &sites)
    : type_(type), matrix_(), coupling_(1.0),
      coupling_name_("HYDRA_COUPLING_NAMELESS"), sites_(sites) {}

Bond::Bond(std::string type, complex coupling, int site)
    : type_(type), matrix_(), coupling_(coupling),
      coupling_name_("HYDRA_COUPLING_NAMELESS"), sites_({site}) {}

Bond::Bond(std::string type, complex coupling, std::vector<int> const &sites)
    : type_(type), matrix_(), coupling_(coupling),
      coupling_name_("HYDRA_COUPLING_NAMELESS"), sites_(sites) {}

Bond::Bond(std::string type, std::string coupling_name, int site)
    : type_(type), matrix_(), coupling_(0.0), coupling_name_(coupling_name),
      sites_({site}) {}

Bond::Bond(std::string type, std::string coupling_name,
           std::vector<int> const &sites)
    : type_(type), matrix_(), coupling_(0.0), coupling_name_(coupling_name),
      sites_(sites) {}

// Constructors with bond matrix
template <typename coeff_t>
Bond::Bond(arma::Mat<coeff_t> const &matrix, int site)
    : type_("HYDRA_TYPE_UNDEFINED"), matrix_(to_cx_mat(matrix)), coupling_(1.0),
      coupling_name_("HYDRA_COUPLING_NAMELESS"), sites_({site}) {}
template Bond::Bond(arma::Mat<double> const &, int);
template Bond::Bond(arma::Mat<complex> const &, int);

template <typename coeff_t>
Bond::Bond(arma::Mat<coeff_t> const &matrix, std::vector<int> const &sites)
    : type_("HYDRA_TYPE_UNDEFINED"), matrix_(to_cx_mat(matrix)), coupling_(1.0),
      coupling_name_("HYDRA_COUPLING_NAMELESS"), sites_(sites) {}
template Bond::Bond(arma::Mat<double> const &, std::vector<int> const &);
template Bond::Bond(arma::Mat<complex> const &, std::vector<int> const &);

template <typename coeff_t>
Bond::Bond(arma::Mat<coeff_t> const &matrix, complex coupling, int site)
    : type_("HYDRA_TYPE_UNDEFINED"), matrix_(to_cx_mat(matrix)),
      coupling_(coupling), coupling_name_("HYDRA_COUPLING_NAMELESS"),
      sites_({site}) {}
template Bond::Bond(arma::Mat<double> const &, complex, int);
template Bond::Bond(arma::Mat<complex> const &, complex, int);

template <typename coeff_t>
Bond::Bond(arma::Mat<coeff_t> const &matrix, complex coupling,
           std::vector<int> const &sites)
    : type_("HYDRA_TYPE_UNDEFINED"), matrix_(to_cx_mat(matrix)),
      coupling_(coupling), coupling_name_("HYDRA_COUPLING_NAMELESS"),
      sites_(sites) {}
template Bond::Bond(arma::Mat<double> const &, complex,
                    std::vector<int> const &);
template Bond::Bond(arma::Mat<complex> const &, complex,
                    std::vector<int> const &);

template <typename coeff_t>
Bond::Bond(arma::Mat<coeff_t> const &matrix, std::string coupling_name,
           int site)
    : type_("HYDRA_TYPE_UNDEFINED"), matrix_(to_cx_mat(matrix)), coupling_(0.0),
      coupling_name_(coupling_name), sites_({site}) {}
template Bond::Bond(arma::Mat<double> const &, std::string, int);
template Bond::Bond(arma::Mat<complex> const &, std::string, int);

template <typename coeff_t>
Bond::Bond(arma::Mat<coeff_t> const &matrix, std::string coupling_name,
           std::vector<int> const &sites)
    : type_("HYDRA_TYPE_UNDEFINED"), matrix_(to_cx_mat(matrix)), coupling_(0.0),
      coupling_name_(coupling_name), sites_(sites) {}
template Bond::Bond(arma::Mat<double> const &, std::string coupling_name,
                    std::vector<int> const &);
template Bond::Bond(arma::Mat<complex> const &, std::string coupling_name,
                    std::vector<int> const &);

bool Bond::type_defined() const { return (type_ != "HYDRA_TYPE_UNDEFINED"); }
bool Bond::matrix_defined() const { return !type_defined(); }
bool Bond::coupling_defined() const {
  return (coupling_name_ == "HYDRA_COUPLING_NAMELESS");
}
bool Bond::coupling_named() const { return !coupling_defined(); }
bool Bond::sites_disjoint() const {
  auto set = std::set<int>(sites_.begin(), sites_.end());
  return set.size() == sites_.size();
}

std::string Bond::type() const {
  if (!type_defined()) {
    Log.err("Error: cannot get type of bond. Type is undefined for this bond.");
  }
  return type_;
}

arma::cx_mat Bond::matrix() const {
  if (!matrix_defined()) {
    Log.err(
        "Error: cannot get matrix of bond. Matrix is undefined for this bond.");
  }
  return matrix_;
}

arma::mat Bond::matrix_real() const {
  if (!matrix_defined()) {
    Log.err("Error: cannot get matrix (real part) of bond. Matrix is undefined "
            "for this bond.");
  }
  return arma::real(matrix_);
}

template <typename coeff_t> coeff_t Bond::coupling(double precision) const {

  if (!coupling_defined()) {
    Log.err("Error: only coupling_name is defined for this bond: {}",
            coupling_name_);
  }

  if constexpr (hydra::is_real<coeff_t>()) {
    if (std::abs(imag(coupling_)) > precision) {
      Log.err("Error: cannot return real coupling for bond. Imaginary part "
              "non-negligible.");
      return 0.;
    }
    return real(coupling_);
  } else {
    (void)precision;
    return coupling_;
  }
}
template double Bond::coupling<double>(double) const;
template complex Bond::coupling<complex>(double) const;

std::string Bond::coupling_name() const {
  if (!coupling_named()) {
    Log.err("Error: no coupling_name is defined for this bond!");
  }
  return coupling_name_;
}
std::vector<int> Bond::sites() const { return sites_; }

int Bond::site(int j) const { return sites_.at(j); }
int Bond::size() const { return (int)sites_.size(); }
int Bond::operator[](int j) const { return site(j); }

bool Bond::is_complex(double precision) const {
  if (type_defined()) {

    // Inherently complex interactions
    if (std::find(complex_bond_types.begin(), complex_bond_types.end(),
                  type_) != complex_bond_types.end()) {
      return true;
    }

    // if coupling is defined, check whether it has imag. part
    if (coupling_defined()) {
      if (std::abs(imag(coupling_)) > precision) {
        return true;
      } else {
        return false;
      }
    } else {
      Log.err(
          "Error: cannot determine if bond is complex if coupling is named!");
    }

    // By default we assume bonds are real
    return false;

    // matrix defined
  } else {
    bool matrix_real = (arma::norm(arma::imag(matrix_)) < precision);

    // if coupling is defined, check whether it has imag. part
    if (coupling_defined()) {
      bool coupling_real = (std::abs(imag(coupling_)) < precision);
      return !(matrix_real && coupling_real);
    } else {
      Log.err(
          "Error: cannot determine if bond is complex if coupling is named!");
    }

    // Dummy, should not happen
    return matrix_real;
  }
}

bool Bond::is_real(double precision) const { return !is_complex(precision); }
bool Bond::operator==(const Bond &rhs) const {
  return (type_ == rhs.type_) &&
         arma::approx_equal(matrix_, rhs.matrix_, "both", 1e-12, 1e-12) &&
         (close(coupling_, rhs.coupling_)) &&
         (coupling_name_ == rhs.coupling_name_) && (sites_ == rhs.sites_);
}

std::vector<int> common_sites(Bond b1, Bond b2) {
  std::vector<int> s1 = b1.sites();
  std::vector<int> s2 = b2.sites();
  std::vector<int> s12;
  std::sort(s1.begin(), s1.end());
  std::sort(s2.begin(), s2.end());
  std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::back_inserter(s12));
  return s12;
}

std::ostream &operator<<(std::ostream &out, const Bond &bond) {

  if (bond.coupling_named()) {
    out << bond.type() << " " << bond.coupling_name() << " ";
  } else {
    out << bond.type() << " " << bond.coupling() << " ";
  }
  for (auto site : bond.sites())
    out << site << " ";
  return out;
}

} // namespace hydra
