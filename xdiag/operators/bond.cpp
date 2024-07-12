#include "bond.hpp"

#include <algorithm>
#include <set>

#include <xdiag/common.hpp>
#include <xdiag/utils/close.hpp>

namespace xdiag {

Bond::Bond(std::string type, coupling_t coupling,
           std::vector<int64_t> const &sites) try
    : type_(type), coupling_(coupling), sites_(sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Bond::Bond(std::string type, coupling_t coupling, int64_t site) try
    : type_(type), coupling_(coupling), sites_({site}) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::string Bond::type() const { return type_; }

template <typename T> bool Bond::coupling_is() const {
  return std::holds_alternative<T>(coupling_);
};
template bool Bond::coupling_is<std::string>() const;
template bool Bond::coupling_is<double>() const;
template bool Bond::coupling_is<complex>() const;
template bool Bond::coupling_is<arma::mat>() const;
template bool Bond::coupling_is<arma::cx_mat>() const;

template <typename T> T Bond::coupling() const try {
  return std::get<T>(coupling_);
} catch (std::bad_variant_access const &ex) {
  XDIAG_THROW("Invalid type given when retrieving coupling");
  return T();
}
template std::string Bond::coupling<std::string>() const;
template double Bond::coupling<double>() const;
template arma::mat Bond::coupling<arma::mat>() const;

template <> complex Bond::coupling<complex>() const try {
  if (coupling_is<double>()) {
    return complex(coupling<double>(), 0.0);
  } else {
    return std::get<complex>(coupling_);
  }
} catch (std::bad_variant_access const &ex) {
  XDIAG_THROW("Invalid type given when retrieving coupling");
  return complex();
}

template <> arma::cx_mat Bond::coupling<arma::cx_mat>() const try {
  if (coupling_is<arma::mat>()) {
    arma::mat mat_r = coupling<arma::mat>();
    arma::mat mat_i = arma::mat(mat_r.n_rows, mat_r.n_cols, arma::fill::zeros);
    return arma::cx_mat(mat_r, mat_i);
  } else {
    std::get<arma::cx_mat>(coupling_);
  }
} catch (std::bad_variant_access const &ex) {
  XDIAG_THROW("Invalid type given when retrieving coupling");
  return arma::cx_mat();
}
coupling_t Bond::coupling() const { return coupling_; }
std::string Bond::coupling_type_string() const try {
  if (std::holds_alternative<std::string>(coupling_)) {
    return "std::string";
  } else if (std::holds_alternative<double>(coupling_)) {
    return "double";
  } else if (std::holds_alternative<complex>(coupling_)) {
    return "complex";
  } else if (std::holds_alternative<arma::mat>(coupling_)) {
    return "arma::mat";
  } else if (std::holds_alternative<arma::cx_mat>(coupling_)) {
    return "arma::cx_mat";
  } else {
    XDIAG_THROW("Invalid coupling type found");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return "";
}

int64_t Bond::size() const { return sites_.size(); }
int64_t Bond::operator[](int64_t idx) const try {
  return sites_.at(idx);
} catch (std::out_of_range const &exc) {
  XDIAG_THROW("Site index out of range for Bond");
  return 0;
}
std::vector<int64_t> const &Bond::sites() const { return sites_; }

bool Bond::isexplicit() const try {
  return !coupling_is<std::string>();
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return false;
}

bool Bond::isreal() const try {
  if (coupling_is<std::string>()) {
    XDIAG_THROW("Cannot determine whether coupling is real, since coupling is "
                "set to a std::string value");
  } else {
    return coupling_is<double>() || coupling_is<arma::mat>();
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return false;
}

bool Bond::ismatrix() const {
  return coupling_is<arma::mat>() || coupling_is<arma::cx_mat>();
}

bool Bond::operator==(const Bond &rhs) const {
  return (type_ == rhs.type_) && (coupling_ == rhs.coupling_) &&
         (sites_ == rhs.sites_);
}
bool Bond::operator!=(const Bond &rhs) const { return !operator==(rhs); }

bool sites_disjoint(Bond const &bond) {
  auto const &sites = bond.sites();
  auto set = std::set<int64_t>(sites.begin(), sites.end());
  return set.size() == sites.size();
}

std::vector<int64_t> common_sites(Bond const &b1, Bond const &b2) {
  std::vector<int64_t> s1 = b1.sites();
  std::vector<int64_t> s2 = b2.sites();
  std::vector<int64_t> s12;
  std::sort(s1.begin(), s1.end());
  std::sort(s2.begin(), s2.end());
  std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::back_inserter(s12));
  return s12;
}

std::ostream &operator<<(std::ostream &out, const Bond &bond) try {

  if (bond.coupling_is<std::string>()) {
    out << bond.type() << " " << bond.coupling<std::string>() << " ";
  } else if (bond.coupling_is<double>()) {
    out << bond.type() << " " << bond.coupling<double>() << " ";
  } else if (bond.coupling_is<complex>()) {
    out << bond.type() << " " << bond.coupling<complex>() << " ";
  } else if (bond.coupling_is<arma::mat>()) {
    out << bond.type() << " " << bond.coupling<arma::mat>() << " ";
  } else if (bond.coupling_is<arma::cx_mat>()) {
    out << bond.type() << " " << bond.coupling<arma::cx_mat>() << " ";
  } else {
    XDIAG_THROW("Unable to print Bond");
  }

  for (int64_t site : bond.sites()) {
    out << site << " ";
  }
  return out;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return out;
}


} // namespace xdiag
