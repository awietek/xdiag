#pragma once

#include <hydra/symmetries/representation.h>

namespace hydra::utils {

template <class bit_t, class GroupAction>
double symmetrized_norm_spinhalf(bit_t spins, GroupAction &&group_action,
                                 Representation const &irrep) {
  assert(group_action.n_symmetries() == irrep.size());
  complex amplitude = 0.0;
  for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tspins = group_action.apply(sym, spins);
    if (tspins == spins) {
      amplitude += irrep.character(sym);
    }
  }
  return std::sqrt(std::abs(amplitude));
}

template <class bit_t, class GroupAction>
double symmetrized_norm_electron(bit_t ups, bit_t dns,
                                 GroupAction &&group_action,
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

} // namespace hydra::utils
