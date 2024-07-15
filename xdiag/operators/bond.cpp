#include "bond.hpp"

#include <algorithm>
#include <set>

#include <xdiag/common.hpp>
#include <xdiag/utils/close.hpp>

namespace xdiag {

Bond::Bond(std::string type, Coupling coupling,
           std::vector<int64_t> const &sites) try
    : type_(type), coupling_(coupling), sites_(sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Bond::Bond(std::string type, Coupling coupling, int64_t site) try
    : type_(type), coupling_(coupling), sites_({site}) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Bond::Bond(std::string type, const char *coupling,
           std::vector<int64_t> const &sites) try
    : Bond(type, Coupling(std::string(coupling)), sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Bond::Bond(std::string type, const char *coupling, int64_t site) try
    : Bond(type, coupling, std::vector<int64_t>({site})) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Bond::Bond(std::string type, std::string coupling,
           std::vector<int64_t> const &sites) try
    : Bond(type, Coupling(coupling), sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Bond::Bond(std::string type, std::string coupling, int64_t site) try
    : Bond(type, coupling, std::vector<int64_t>({site})) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Bond::Bond(std::string type, double coupling,
           std::vector<int64_t> const &sites) try
    : Bond(type, Coupling(coupling), sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Bond::Bond(std::string type, double coupling, int64_t site) try
    : Bond(type, coupling, std::vector<int64_t>({site})) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Bond::Bond(std::string type, complex coupling,
           std::vector<int64_t> const &sites) try
    : Bond(type, Coupling(coupling), sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Bond::Bond(std::string type, complex coupling, int64_t site) try
    : Bond(type, coupling, std::vector<int64_t>({site})) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Bond::Bond(std::string type, arma::mat const &coupling,
           std::vector<int64_t> const &sites) try
    : Bond(type, Coupling(coupling), sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Bond::Bond(std::string type, arma::mat const &coupling, int64_t site) try
    : Bond(type, coupling, std::vector<int64_t>({site})) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Bond::Bond(std::string type, arma::cx_mat const &coupling,
           std::vector<int64_t> const &sites) try
    : Bond(type, Coupling(coupling), sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Bond::Bond(std::string type, arma::cx_mat const &coupling, int64_t site) try
    : Bond(type, coupling, std::vector<int64_t>({site})) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::string Bond::type() const { return type_; }

Coupling const &Bond::coupling() const { return coupling_; }

int64_t Bond::size() const { return sites_.size(); }
int64_t Bond::operator[](int64_t idx) const try {
  return sites_.at(idx);
} catch (std::out_of_range const &exc) {
  XDIAG_THROW("Site index out of range for Bond");
}
std::vector<int64_t> const &Bond::sites() const { return sites_; }

bool Bond::isexplicit() const try {
  return coupling_.isexplicit();
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool Bond::isreal() const try {
  if (coupling_.is<std::string>()) {
    XDIAG_THROW("Cannot determine whether coupling is real, since coupling is "
                "set to a std::string value");
  } else if (type_ == "SCALARCHIRALITY") {
    return false;
  } else {
    return coupling_.isreal();
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool Bond::ismatrix() const { return coupling_.ismatrix(); }

bool Bond::operator==(Bond const &rhs) const try {
  return (type_ == rhs.type_) && (coupling_ == rhs.coupling_) &&
         (sites_ == rhs.sites_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool Bond::operator!=(Bond const &rhs) const { return !operator==(rhs); }

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

  out << bond.type() << " " << bond.coupling() << " ";
  for (int64_t site : bond.sites()) {
    out << site << " ";
  }

  return out;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
