#pragma once

#include <string>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra::utils {

bool coupling_is_zero(Bond const &bond, Couplings const &couplings);
bool coupling_is_non_zero(Bond const &bond, Couplings const &couplings);

void check_operator_real(BondList const &bonds, Couplings const &cpls,
                         std::string errmsg);

void check_operator_real(BondList const &bonds, Couplings const &cpls,
                         Representation const &irrep_in,
                         Representation const &irrep_out, std::string errmsg);

template <typename coeff_t>
void check_operator_works_with(BondList const &bonds, Couplings const &cpls,
                               std::string errmsg);
template <typename coeff_t>
void check_operator_works_with(BondList const &bonds, Couplings const &cpls,
                               Representation const &irrep_in,
                               Representation const &irrep_out,
                               std::string errmsg);

bool sites_disjoint(std::vector<int> const &sites);
bool sites_disjoint(Bond const &bond);

void check_sites_disjoint(std::vector<int> const &sites);
void check_sites_disjoint(Bond const &bond);

  
BondList clean_bondlist(BondList const &bonds, Couplings const &cpls,
                        std::vector<std::string> desired_bond_types);
BondList clean_bondlist(BondList const &bonds, Couplings const &cpls,
                        std::vector<std::string> desired_bond_types,
                        int allowed_size);
BondList non_zero_bondlist(BondList const &bonds, Couplings const &cpls);

template <class coeff_t>
coeff_t get_coupling(Couplings const &couplings, std::string cpl);

template <class coeff_t>
std::pair<coeff_t, coeff_t> get_coupling_and_conj(Couplings const &couplings,
                                                  std::string cpl);

} // namespace hydra::utils
