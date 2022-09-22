#pragma once

#include <string>

#include <hydra/operators/bondlist.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra::utils {

void check_nup_spinhalf(int n_sites, int nup, std::string block_name);
void check_nup_ndn_tj(int n_sites, int nup, int ndn, std::string block_name);
void check_nup_ndn_electron(int n_sites, int nup, int ndn,
                            std::string block_name);

void check_n_sites(int n_sites, PermutationGroup const &permutation_group);

template <class coeff_t>
std::vector<coeff_t> characters(Representation const &irrep);

} // namespace hydra::utils
