#pragma once

#include <string>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/representation.h>

namespace hydra::utils {

bool coupling_is_zero(Bond const &bond, Couplings const &couplings);
bool coupling_is_non_zero(Bond const &bond, Couplings const &couplings);
void check_nup_spinhalf(int n_sites, int nup, std::string block_name);
void check_nup_ndn_tj(int n_sites, int nup, int ndn, std::string block_name);
void check_nup_ndn_electron(int n_sites, int nup, int ndn,
                            std::string block_name);

void check_operator_real(BondList const &bonds, Couplings const &cpls,
                         std::string errmsg);
void check_symmetric_operator_real(BondList const &bonds, Couplings const &cpls,
                                   Representation const &irrep_in,
                                   Representation const &irrep_out,
                                   std::string errmsg);

BondList clean_bondlist(BondList const &bonds, Couplings const &cpls,
                        std::vector<std::string> desired_bond_types,
                        int allowed_size);
void warn_if_complex(complex x, std::string msg = "");

} // namespace hydra::utils
