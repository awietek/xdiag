#pragma once

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <hydra/symmetries/representation.h>

namespace hydra::utils {

bool coupling_is_zero(Bond const &bond, Couplings const &couplings);

template <class bit_t, class GroupAction>
double compute_norm(bit_t ups, bit_t dns, GroupAction &&group_action,
                    Representation const &irrep) {
  assert(group_action.n_symmetries() == irrep.size());
  complex amplitude = 0.0;
  for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tups = group_action.apply(sym, ups);
    if (tups == ups) {
      bit_t tdns = group_action.apply(sym, dns);
      double fermi_sign_ups = group_action.fermi_sign(sym, ups);
      if (tdns == dns) {
        double fermi_sign_dns = group_action.fermi_sign(sym, dns);
        amplitude += fermi_sign_ups * fermi_sign_dns * irrep.character(sym);
      }
    }
  }
  return std::sqrt(std::abs(amplitude));
}

void check_nup_ndn_tj(int n_sites, int nup, int ndn);
void check_nup_ndn_electron(int n_sites, int nup, int ndn);


} // namespace hydra::utils
