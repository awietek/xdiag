#pragma once

#include <string>

#include <xdiag/operators/bondlist.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::utils {

void check_nup_spinhalf(int nsites, int nup, std::string block_name);
void check_nup_ndn_tj(int nsites, int nup, int ndn, std::string block_name);
void check_nup_ndn_electron(int nsites, int nup, int ndn,
                            std::string block_name);

void check_nsites(int nsites, PermutationGroup const &permutation_group);

template <class coeff_t>
std::vector<coeff_t> characters(Representation const &irrep);

} // namespace xdiag::utils
