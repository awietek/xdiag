#include "bondlist.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

#include <hydra/utils/print_macro.h>

namespace hydra {

BondList::BondList(std::vector<Bond> const &bonds) : bonds_(bonds) {}
BondList::BondList(io::FileTomlHandler &&hdl) : BondList(hdl.as<BondList>()) {}

void BondList::operator<<(Bond const &bond) { bonds_.push_back(bond); }

bool BondList::coupling_defined(std::string name) const {
  return couplings_.count(name);
}
void BondList::set_coupling(std::string name, complex cpl) {
  couplings_[name] = cpl;
}

template <typename coeff_t>
coeff_t BondList::coupling(std::string name, double precision) const {
  if (couplings_.count(name)) {
    complex cpl = couplings_.at(name);
    if constexpr (hydra::isreal<coeff_t>()) {
      if (std::abs(imag(cpl)) > precision) {
        Log.err("Error: cannot return real coupling for bond. Imaginary part "
                "non-negligible.");
        return 0.;
      }
      return real(cpl);
    } else {
      (void)precision;
      return cpl;
    }
  } else {
    Log.err("Error: undefined coupling in BondList: {}", name);
    return 0.;
  }
}
template double BondList::coupling<double>(std::string, double) const;
template complex BondList::coupling<complex>(std::string, double) const;

std::map<std::string, complex> const &BondList::couplings() const {
  return couplings_;
}

bool BondList::matrix_defined(std::string name) const {
  return matrices_.count(name);
}
void BondList::set_matrix(std::string name, arma::cx_mat mat) {
  matrices_[name] = mat;
}
void BondList::set_matrix(std::string name, arma::mat mat) {
  matrices_[name] = to_cx_mat(mat);
}
arma::cx_mat BondList::matrix(std::string name) const {
  if (matrices_.count(name)) {
    return matrices_.at(name);
  } else {
    Log.err("Error: undefined matrix in BondList: {}", name);
  }
  return matrices_.at(name);
}
std::map<std::string, arma::cx_mat> const &BondList::matrices() const {
  return matrices_;
}

BondListHandler BondList::operator[](std::string name) {
  return BondListHandler(name, couplings_, matrices_);
}

int64_t BondList::size() const { return (int64_t)bonds_.size(); }
void BondList::clear() {
  bonds_.clear();
  couplings_.clear();
  matrices_.clear();
}

int64_t BondList::n_sites() const {
  int64_t n_sites = 0;
  for (Bond bond : bonds_) {
    for (int64_t site : bond.sites()) {
      n_sites = std::max(n_sites, site + 1);
    }
  }
  return n_sites;
}

BondList BondList::bonds_of_type(std::string type) const {
  BondList bonds_return;
  for (Bond bond : bonds_) {
    if ((bond.type_defined()) && (bond.type() == type)) {
      if (bond.coupling_named()) {
        std::string name = bond.coupling_name();
        bonds_return << Bond(type, name, bond.sites());

        if (coupling_defined(name)) {
          complex cpl = coupling(name);
          bonds_return[name] = cpl;
        }
      } else {
        bonds_return << Bond(type, bond.coupling(), bond.sites());
      }
    }
  }
  return bonds_return;
}

BondList BondList::bonds_with_matrix() const {
  std::vector<Bond> bonds_return;
  for (Bond bond : bonds_) {
    if (bond.matrix_defined()) {
      bonds_return.push_back(bond);
    }
  }
  return BondList(bonds_return);
}

bool BondList::iscomplex(double precision) const {
  for (Bond bond : bonds_) {

    if (bond.matrix_defined()) {
      auto mat = bond.matrix();
      bool matrix_real = std::abs(arma::norm(arma::imag(mat))) < precision;

      bool coupling_real = false;
      if (bond.coupling_defined()) {
        coupling_real = std::abs(imag(bond.coupling())) < precision;
      } else {
        std::string coupling_name = bond.coupling_name();

        // coupling name defined in couplings
        if (couplings_.count(coupling_name)) {
          complex cpl = couplings_.at(coupling_name);
          coupling_real = std::abs(imag(cpl)) < precision;
        } else {
          Log.err("Error: cannot determine whether Bond in BondList is "
                  "complex/real. Its coupling coefficient is not uniquely "
                  "determined.");
        }
      }
      if ((!matrix_real) || (!coupling_real)) {
        return true;
      }

    } else { // bond.type_defined()

      // Inherently complex interactions
      if (std::find(complex_bond_types.begin(), complex_bond_types.end(),
                    bond.type()) != complex_bond_types.end()) {
        return true;
      }

      // if coupling is defined, check whether it has imag. part
      bool coupling_real = false;
      if (bond.coupling_defined()) {
        coupling_real = std::abs(imag(bond.coupling())) < precision;
      } else {
        std::string coupling_name = bond.coupling_name();

        // coupling name defined in couplings
        if (couplings_.count(coupling_name)) {
          complex cpl = couplings_.at(coupling_name);
          coupling_real = (std::abs(imag(cpl)) < precision);
        } else {
          Log.err("Error: cannot determine whether Bond in BondList is "
                  "complex/real. Its coupling coefficient is not uniquely "
                  "determined.");
        }
      }
      if (!coupling_real) {
        return true;
      }
    }
  }
  return false;
}

BondList::iterator_t BondList::begin() { return bonds_.begin(); }
BondList::iterator_t BondList::end() { return bonds_.end(); }
BondList::const_iterator_t BondList::begin() const { return bonds_.begin(); }
BondList::const_iterator_t BondList::end() const { return bonds_.end(); }
BondList::const_iterator_t BondList::cbegin() const { return bonds_.cbegin(); }
BondList::const_iterator_t BondList::cend() const { return bonds_.cend(); }

bool BondList::isreal(double precision) const { return !iscomplex(precision); }
bool BondList::ishermitian(double precision) const {
  (void) precision;
  return true; // implement this
}

bool BondList::operator==(BondList const &other) const {
  std::vector<std::string> matrices_keys;
  for (const auto &[key, val] : matrices_) {
    (void)val;
    matrices_keys.push_back(key);
  }

  std::vector<std::string> other_matrices_keys;
  for (const auto &[key, _] : other.matrices_) {
    (void)_;
    other_matrices_keys.push_back(key);
  }

  if (matrices_keys != other_matrices_keys) {
    return false;
  } else {
    for (auto &&key : matrices_keys) {
      if (!approx_equal(matrices_.at(key), other.matrices_.at(key), "both",
                        1e-12, 1e-12)) {
        return false;
      }
    }
  }
  return (bonds_ == other.bonds_) && (couplings_ == other.couplings_);
}

bool BondList::operator!=(BondList const &other) const {
  return !operator==(other);
}

BondList BondList::operator+(BondList const &other) {
  auto newbonds = bonds_;
  newbonds.insert(newbonds.end(), other.begin(), other.end());
  BondList blnew = BondList(newbonds);

  for (auto it : matrices_) {
    blnew.set_matrix(it.first, it.second);
  }

  for (auto it : couplings_) {
    blnew.set_coupling(it.first, it.second);
  }

  for (auto it : other.matrices()) {
    blnew.set_matrix(it.first, it.second);
  }

  for (auto it : other.couplings()) {
    blnew.set_coupling(it.first, it.second);
  }

  return blnew;
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

} // namespace hydra
