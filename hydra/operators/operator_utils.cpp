#include "operator_utils.h"

#include <hydra/utils/complex.h>
#include <hydra/utils/logger.h>

#include <set>

namespace hydra::utils {

bool coupling_is_zero(Bond const &bond, Couplings const &couplings) {
  std::string coupling = bond.coupling();
  return (!couplings.defined(coupling)) ||
         lila::close(couplings[coupling], (complex)0.);
}

bool coupling_is_non_zero(Bond const &bond, Couplings const &couplings) {
  return !coupling_is_zero(bond, couplings);
}
void check_n_sites(int n_sites, PermutationGroup const &permutation_group) {
  if (n_sites != permutation_group.n_sites()) {
    Log.err(
        "Error creating block: n_sites disagrees with permutation_group");
  }
}

void check_operator_real(BondList const &bonds, Couplings const &cpls,
                         std::string errmsg) {
  if (is_complex(bonds))
    Log.err("Error: cannot {} from complex bonds!", errmsg);
  if (is_complex(cpls))
    Log.err("Error: cannot {} real matrix from complex couplings!",
                  errmsg);
}

void check_operator_real(BondList const &bonds, Couplings const &cpls,
                         Representation const &irrep_in,
                         Representation const &irrep_out, std::string errmsg) {
  check_operator_real(bonds, cpls, errmsg);
  if (is_complex(irrep_in) || is_complex(irrep_out))
    Log.err("Error: cannot {} real matrix from complex representation!",
                  errmsg);
}

template <typename coeff_t>
void check_operator_works_with(BondList const &bonds, Couplings const &cpls,
                               std::string errmsg) {
  if constexpr (is_real<coeff_t>()) {
    check_operator_real(bonds, cpls, errmsg);
  }
}
template void check_operator_works_with<double>(BondList const &,
                                                Couplings const &, std::string);
template void check_operator_works_with<complex>(BondList const &,
                                                 Couplings const &,
                                                 std::string);
template <typename coeff_t>
void check_operator_works_with(BondList const &bonds, Couplings const &cpls,
                               Representation const &irrep_in,
                               Representation const &irrep_out,
                               std::string errmsg) {
  if constexpr (is_real<coeff_t>()) {
    check_operator_real(bonds, cpls, irrep_in, irrep_out, errmsg);
  }
}
template void check_operator_works_with<double>(BondList const &,
                                                Couplings const &,
                                                Representation const &,
                                                Representation const &,
                                                std::string);
template void check_operator_works_with<complex>(BondList const &,
                                                 Couplings const &,
                                                 Representation const &,
                                                 Representation const &,
                                                 std::string);

void check_sites_disjoint(std::vector<int> const &sites) {
  auto set = std::set<int>(sites.begin(), sites.end());
  if (set.size() != sites.size())
    Log.err("Error: sites are not disjoint");
}

void check_sites_disjoint(Bond const &bond) {
  check_sites_disjoint(bond.sites());
}


BondList clean_bondlist(BondList const &bonds, Couplings const &cpls,
                        std::vector<std::string> desired_bond_types) {
  BondList clean_bonds;

  for (auto type : desired_bond_types) {
    auto bonds_of_type = bonds.bonds_of_type(type);

    for (auto bond : bonds_of_type) {

      if (coupling_is_non_zero(bond, cpls)) {
        clean_bonds << bond;
      }
    }
  }
  return clean_bonds;
}

BondList clean_bondlist(BondList const &bonds, Couplings const &cpls,
                        std::vector<std::string> desired_bond_types,
                        int allowed_size) {
  BondList clean_bonds;

  for (auto type : desired_bond_types) {
    auto bonds_of_type = bonds.bonds_of_type(type);

    for (auto bond : bonds_of_type) {
      if (bond.size() != allowed_size) {
        Log.err("Invalid bond size ({}) found for type {}", bond.size(),
		type);
      }

      if (coupling_is_non_zero(bond, cpls)) {
        clean_bonds << bond;
      }
    }
  }
  return clean_bonds;
}

BondList non_zero_bondlist(BondList const &bonds, Couplings const &cpls) {
  BondList clean_bonds;

  for (auto bond : bonds) {

    if (coupling_is_non_zero(bond, cpls)) {
      clean_bonds << bond;
    }
  }
  return clean_bonds;
}

template <class coeff_t>
coeff_t get_coupling(Couplings const &couplings, std::string cpl) {
  coeff_t val;
  if constexpr (is_complex<coeff_t>()) {
    val = couplings[cpl];
  } else {
    utils::warn_if_complex(couplings[cpl], "imaginary part deprecated");
    val = lila::real(couplings[cpl]);
  }
  return val;
}

template double get_coupling<double>(Couplings const &couplings,
                                     std::string cpl);
template complex get_coupling<complex>(Couplings const &couplings,
                                       std::string cpl);

template <class coeff_t>
std::pair<coeff_t, coeff_t> get_coupling_and_conj(Couplings const &couplings,
                                                  std::string cpl) {
  coeff_t val;
  coeff_t val_conj;
  if constexpr (is_complex<coeff_t>()) {
    val = couplings[cpl];
    val_conj = lila::conj(val);

  } else {
    utils::warn_if_complex(couplings[cpl], "imaginary part deprecated");
    val = lila::real(couplings[cpl]);
    val_conj = val;
  }
  return {val, val_conj};
}

template std::pair<double, double>
get_coupling_and_conj<double>(Couplings const &couplings, std::string cpl);
template std::pair<complex, complex>
get_coupling_and_conj<complex>(Couplings const &couplings, std::string cpl);

} // namespace hydra::utils
