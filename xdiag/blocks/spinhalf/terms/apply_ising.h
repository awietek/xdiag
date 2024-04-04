#pragma once

#include <xdiag/bits/bitops.h>
#include <xdiag/common.h>

#include <xdiag/blocks/spinhalf/terms/apply_term_diag.h>

namespace xdiag::spinhalf {

// Ising term: J S^z_i S^z_j

template <typename bit_t, typename coeff_t, bool symmetric, class BasisIn,
          class BasisOut, class Fill>
void apply_ising(Bond const &bond, BasisIn &&basis_in,
                 BasisOut &&basis_out, Fill &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined() && (bond.type() == "ISING"));
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  int64_t s1 = bond[0];
  int64_t s2 = bond[1];
  bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

  coeff_t J = bond.coupling<coeff_t>();

  coeff_t val_same = J / 4.0;
  coeff_t val_diff = -J / 4.0;

  if (basis_in == basis_out) {

    auto term_coeff = [&mask, &val_same, &val_diff](bit_t spins) -> coeff_t {
      if (bits::popcnt(spins & mask) & 1) {
        return val_diff;
      } else {
        return val_same;
      }
    };
    spinhalf::apply_term_diag<bit_t, coeff_t>(basis_in, term_coeff, fill);

  } else {
    auto non_zero_term = [](bit_t spins) -> bool {
      return true;
      (void)spins;
    };
    auto term_action = [&mask, &val_same,
                        &val_diff](bit_t spins) -> std::pair<bit_t, coeff_t> {
      if (bits::popcnt(spins & mask) & 1) {
        return {spins, val_diff};
      } else {
        return {spins, val_same};
      }
    };
    if constexpr (symmetric) {
      spinhalf::apply_term_offdiag_sym<bit_t, coeff_t>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    } else {
      spinhalf::apply_term_offdiag_no_sym<bit_t, coeff_t>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    }
  }
}

} // namespace xdiag::spinhalf
