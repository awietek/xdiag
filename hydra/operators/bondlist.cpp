#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

#include "bondlist.h"

namespace hydra {

namespace detail {

std::vector<TypeCoupling> get_types_couplings(std::vector<Bond> const &bonds) {
  std::vector<TypeCoupling> types_couplings;

  for (Bond bond : bonds) {
    TypeCoupling tc = type_coupling(bond);
    if (std::find(types_couplings.begin(), types_couplings.end(), tc) ==
        types_couplings.end())
      types_couplings.push_back(tc);
  }
  return types_couplings;
}

std::vector<std::string> get_types(std::vector<Bond> const &bonds) {
  std::vector<std::string> types;
  for (Bond bond : bonds)
    if (std::find(types.begin(), types.end(), bond.type()) == types.end())
      types.push_back(bond.type());
  return types;
}

std::vector<std::string> get_couplings(std::vector<Bond> const &bonds) {
  std::vector<std::string> couplings;
  for (Bond bond : bonds)
    if (std::find(couplings.begin(), couplings.end(), bond.coupling()) ==
        couplings.end())
      couplings.push_back(bond.coupling());
  return couplings;
}

std::vector<Bond> get_bonds_of_type(std::vector<Bond> const &bonds,
                                    std::string type) {
  std::vector<Bond> bonds_return;
  for (Bond bond : bonds)
    if (bond.type() == type)
      bonds_return.push_back(bond);
  return bonds_return;
}

std::vector<Bond> get_bonds_of_coupling(std::vector<Bond> const &bonds,
                                        std::string coupling) {
  std::vector<Bond> bonds_return;
  for (Bond bond : bonds)
    if (bond.coupling() == coupling)
      bonds_return.push_back(bond);
  return bonds_return;
}

std::vector<Bond>
get_bonds_of_type_coupling(std::vector<Bond> const &bonds,
                           TypeCoupling const &type_coupling) {
  std::vector<Bond> bonds_return;
  for (Bond bond : bonds)
    if ((bond.type() == type_coupling.type()) &&
        (bond.coupling() == type_coupling.coupling()))
      bonds_return.push_back(bond);
  return bonds_return;
}

} // namespace detail

BondList::BondList(std::vector<Bond> const &bonds) : bonds_(bonds) {}

int BondList::n_sites() const {
  int n_sites = 0;
  for (Bond bond : bonds_)
    for (int site : bond.sites())
      n_sites = std::max(n_sites, site + 1);

  return n_sites;
}

void BondList::operator<<(Bond const &bond) { bonds_.push_back(bond); }

std::vector<std::string> BondList::types() const {
  return detail::get_types(bonds_);
}
std::vector<std::string> BondList::couplings() const {
  return detail::get_couplings(bonds_);
}
std::vector<TypeCoupling> BondList::types_couplings() const {
  return detail::get_types_couplings(bonds_);
}

BondList BondList::bonds_of_type(std::string type) const {
  return BondList(detail::get_bonds_of_type(bonds_, type));
}

BondList BondList::bonds_of_coupling(std::string coupling) const {
  return BondList(detail::get_bonds_of_coupling(bonds_, coupling));
}

BondList BondList::bonds_of_type_coupling(std::string type,
                                          std::string coupling) const {
  return BondList(
      detail::get_bonds_of_type_coupling(bonds_, TypeCoupling(type, coupling)));
}

BondList operator+(BondList const &bl1, BondList const &bl2) {
  auto newbonds = bl1.bonds_;
  newbonds.insert(newbonds.end(), bl2.begin(), bl2.end());
  return BondList(newbonds);
}

BondList read_bondlist(std::string filename) {
  std::vector<Bond> bonds;

  // Open file and handle error
  std::ifstream File(filename.c_str());
  if (File.fail()) {
    std::cerr << "Error in read_bondlist: "
              << "Could not open file with filename [" << filename
              << "] given. Abort." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Advance to interaction lines
  std::string tobeparsed;
  getline(File, tobeparsed);
  while ((tobeparsed.find("[Interactions]") == std::string::npos) &&
         (tobeparsed.find("[interactions]") == std::string::npos))
    getline(File, tobeparsed);

  // read lines until '[' is found or else until EOF
  while (std::getline(File, tobeparsed)) {
    if ((tobeparsed.find('[') != std::string::npos))
      break;

    std::string type, coupling;
    std::vector<int> sites;
    std::stringstream stream(tobeparsed);
    stream >> type;
    stream >> coupling;

    // Just parse sites
    int n;
    while (stream >> n)
      sites.push_back(n);
    bonds.push_back(Bond(type, coupling, sites));

    if (!File.good())
      break;
  }

  return BondList(bonds);
}

bool is_complex(BondList const &bonds) {
  return std::find_if(bonds.begin(), bonds.end(), [](Bond const &bond) {
           return is_complex(bond);
         }) != bonds.end();
}

bool is_real(BondList const &bonds) { return !is_complex(bonds); }

} // namespace hydra
