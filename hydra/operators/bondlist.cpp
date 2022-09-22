#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

#include "bondlist.h"

namespace hydra {

BondList::BondList(std::vector<Bond> const &bonds) : bonds_(bonds) {}

void BondList::operator<<(Bond const &bond) { bonds_.push_back(bond); }

bool BondList::coupling_defined(std::string name) {
  return couplings_.count(name);
}
void BondList::set_coupling(std::string name, complex cpl) {
  couplings_[name] = cpl;
}
complex BondList::get_coupling(std::string name) const {
  if (couplings_.count(name)) {
    return couplings_.at(name);
  } else {
    Log.err("Error: undefined coupling in BondList: {}", name);
    return 0.;
  }
}

bool BondList::matrix_defined(std::string name) {
  return matrices_.count(name);
}
void BondList::set_matrix(std::string name, arma::cx_mat mat) {
  matrices_[name] = mat;
}
void BondList::set_matrix(std::string name, arma::mat mat) {
  matrices_[name] = to_complex(mat);
}
arma::cx_mat BondList::get_matrix(std::string name) {
  if (couplings_.count(name)) {
    return matrices_.at(name);
  } else {
    Log.err("Error: undefined coupling in BondList: {}", name);
    return 0.;
  }
  return matrices_.at(name);
}

int BondList::size() const { return (int)bonds_.size(); }
int BondList::clear() {
  bonds_.clear();
  couplings_.clear();
  matrices_.clear();
}

int BondList::n_sites() const {
  int n_sites = 0;
  for (Bond bond : bonds_) {
    for (int site : bond.sites()) {
      n_sites = std::max(n_sites, site + 1);
    }
  }
  return n_sites;
}

BondList BondList::bonds_of_type(std::string type) const {
  std::vector<Bond> bonds_return;
  for (Bond bond : bonds) {
    if ((bond.type_defined()) && (bond.type() == type)) {
      bonds_return.push_back(bond);
    }
  }
  return bonds_return;
}

BondList BondList::bonds_with_matrix() const {
  std::vector<Bond> bonds_return;
  for (Bond bond : bonds) {
    if (bond.matrix_defined()) {
      bonds_return.push_back(bond);
    }
  }
  return bonds_return;
}

bool BondList::is_complex(double precision) const {
  for (Bond bond : bonds) {

    if (bond.matrix_defined()) {
      auto mat = bond.matrix();
      bool matrix_real = std::abs(arma::norm(arma::imag(mat))) < precision;

      bool coupling_real = false;
      if (bond.coupling_defined()) {
        coupling_real = std::abs(bond.coupling()) < precision;
      } else {
        std::string coupling_name = bond.coupling_name();

        // coupling name defined in couplings
        if (couplings_.count(coupling_name)) {
          complex cpl = couplings_[coupling_name];
          coupling_real = std::abs(cpl) < precision;
        } else {
          Log.err("Error: cannot determine whether Bond in BondList is "
                  "complex/real. Its coupling coefficient is not uniquely "
                  "determined.")
        }
      }
      return matrix_real && coupling_real;
    }

    if (bond.type_defined()) {

      // Inherently complex interactions
      if (std::find(complex_bond_types.begin(), complex_bond_types.end(),
                    type_) != complex_bond_types.end()) {
        return true;
      }

      // if coupling is defined, check whether it has imag. part
      bool coupling_real = false;
      if (bond.coupling_defined()) {
        coupling_real = std::abs(bond.coupling()) < precision;
      } else {
        std::string coupling_name = bond.coupling_name();

        // coupling name defined in couplings
        if (couplings_.count(coupling_name)) {
          complex cpl = couplings_[coupling_name];
          coupling_real = std::abs(cpl) < precision;
        } else {
          Log.err("Error: cannot determine whether Bond in BondList is "
                  "complex/real. Its coupling coefficient is not uniquely "
                  "determined.")
        }
      }
      return coupling_real;
    }
  }
}

bool BondList::is_real(double precision) const {
  return !is_complex(precision);
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

} // namespace hydra
