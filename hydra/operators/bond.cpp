#include "bond.h"

#include <algorithm>

namespace hydra {

Bond::Bond(std::string type, std::string coupling,
           std::vector<int> const &sites)
    : type_(type), coupling_(coupling), sites_(sites), has_parameters_(false),
      parameters_() {}

Bond::Bond(std::string type, std::string coupling,
           std::vector<int> const &sites, Parameters const &parameters)
    : type_(type), coupling_(coupling), sites_(sites), has_parameters_(true),
      parameters_(parameters) {}

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
  out << bond.type() << " " << bond.coupling() << " ";
  for (auto site : bond.sites())
    out << site << " ";
  return out;
}

bool operator==(const Bond &lhs, const Bond &rhs) {
  return (lhs.type() == rhs.type()) && (lhs.coupling() == rhs.coupling()) &&
         (lhs.sites() == rhs.sites());
}

TypeCoupling::TypeCoupling(const std::string &type, const std::string &coupling)
    : type_(type), coupling_(coupling) {}

std::ostream &operator<<(std::ostream &out, const TypeCoupling &tc) {
  out << tc.type() << " " << tc.coupling() << " ";
  return out;
}

bool operator==(const TypeCoupling &tc1, const TypeCoupling &tc2) {
  return (tc1.type() == tc2.type()) && (tc1.coupling() == tc2.coupling());
}

TypeCoupling type_coupling(const Bond &bond) {
  return TypeCoupling(bond.type(), bond.coupling());
}

bool is_complex(Bond const &bond) {
  if (std::find(complex_bond_types.begin(),
                complex_bond_types.end(),
                bond.type()) != complex_bond_types.end()) {
    return true;
  } else {
    return false;
  }
}

bool is_real(Bond const &bond) {
  return !is_complex(bond);
}

} // namespace hydra
