#include "bondlist.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

#include <xdiag/utils/print_macro.hpp>

namespace xdiag {

BondList::BondList(std::vector<Bond> const &bonds) : bonds_(bonds) {}
int64_t BondList::size() const { return bonds_.size(); }
Coupling &BondList::operator[](std::string name) { return couplings_[name]; }
Coupling const &BondList::operator[](std::string name) const {
  return couplings_.at(name);
}

bool BondList::defined(std::string name) const {
  return couplings_.find(name) != couplings_.end();
}

std::vector<std::string> BondList::couplings() const {
  std::vector<std::string> names;
  for (auto const &cpl : couplings_) {
    names.push_back(cpl.first);
  }
  return names;
}

bool BondList::isreal() const try {
  BondList bonds_explicit = make_explicit(*this);
  return std::all_of(bonds_explicit.begin(), bonds_explicit.end(),
                     [](Bond const &b) { return b.isreal(); });
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool BondList::isexplicit() const try {
  return std::all_of(bonds_.begin(), bonds_.end(),
                     [](Bond const &b) { return b.isexplicit(); });
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void BondList::operator+=(Bond const &bond) { bonds_.push_back(bond); }
void BondList::operator+=(BondList const &bonds) try {
  // Add bonds
  for (auto const &bond : bonds.bonds_) {
    bonds_.push_back(bond);
  }

  // Add possible couplings
  for (auto it = bonds.couplings_.begin(); it != bonds.couplings_.end(); it++) {
    std::string name = it->first;
    Coupling cpl = it->second;

    // if name does not exist yet, add coupling
    if (couplings_.find(name) == couplings_.end()) {
      couplings_[name] = cpl;
    }
    // If it exists, throw an error if the coupling disagrees
    else {
      if (couplings_[name] != cpl) {
        XDIAG_THROW(fmt::format("Conflicting values of coupling \"{}\" found "
                                "when trying to add two BondLists.",
                                name));
      }
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

BondList BondList::operator+(Bond const &bond) const {
  BondList new_bonds = *this;
  new_bonds += bond;
  return new_bonds;
}

BondList BondList::operator+(BondList const &bonds) const try {
  BondList new_bonds = *this;
  new_bonds += bonds;
  return new_bonds;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool BondList::operator==(BondList const &other) const {
  if (bonds_ != other.bonds_) {
    return false;
  }
  auto c1s = couplings();
  auto c2s = other.couplings();
  if (c1s != c2s) {
    return false;
  }
  for (std::string c : c1s) {
    if (couplings_.at(c) != other.couplings_.at(c)) {
      return false;
    }
  }
  return true;
}
bool BondList::operator!=(BondList const &other) const {
  return !operator==(other);
}

BondList::iterator_t BondList::begin() { return bonds_.begin(); }
BondList::iterator_t BondList::end() { return bonds_.end(); }
BondList::const_iterator_t BondList::begin() const { return bonds_.begin(); }
BondList::const_iterator_t BondList::end() const { return bonds_.end(); }
BondList::const_iterator_t BondList::cbegin() const { return bonds_.cbegin(); }
BondList::const_iterator_t BondList::cend() const { return bonds_.cend(); }

BondList make_explicit(BondList const &bonds) try {
  BondList explicit_bonds;
  for (auto const &bond : bonds) {
    if (bond.isexplicit()) {
      explicit_bonds += bond;
    } else {
      std::string name = bond.coupling().as<std::string>();
      if (bonds.defined(name)) {
        std::string type = bond.type();
        Coupling cpl = bonds[name];
        if (!cpl.isexplicit()) {
          XDIAG_THROW(fmt::format("Unable to make BondList explicit: coupling "
                                  "\"{}\" is yet a string",
                                  name));
        }
        explicit_bonds += Bond(type, cpl, bond.sites());
      } else {
        XDIAG_THROW(fmt::format(
            "Unable to make BondList explicit: coupling \"{}\" is not defined",
            name));
      }
    }
  }
  return explicit_bonds;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

BondList bonds_of_type(std::string type, BondList const &bonds) {
  BondList new_bonds;
  for (auto const &bond : bonds) {
    if (bond.type() == type) {
      new_bonds += bond;
    }
  }
  return new_bonds;
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
    std::vector<int64_t> sites;
    std::stringstream stream(tobeparsed);
    stream >> type;
    stream >> coupling;

    // Just parse sites
    int64_t n;
    while (stream >> n)
      sites.push_back(n);
    bonds.push_back(Bond(type, coupling, sites));

    if (!File.good())
      break;
  }

  return BondList(bonds);
}

} // namespace xdiag
