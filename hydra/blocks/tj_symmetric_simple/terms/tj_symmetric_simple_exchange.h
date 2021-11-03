#pragma once

#include <lila/utils/logger.h>

#include <hydra/common.h>
#include <hydra/bitops/bitops.h>

#include <hydra/blocks/blocks.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/symmetries/symmetry_utils.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms::tj_symmetric_simple {

template <class bit_t, class coeff_t, class Filler, class GroupAction>
void do_down_flips(bit_t up, idx_t idx_up, bit_t mask, bit_t spacemask,
                   bit_t dnmask, double Jhalf,
                   tJSymmetricSimple<bit_t, GroupAction> const &block,
                   Filler &&fill) {
  using bitops::popcnt;

  auto &group_action = block.group_action();
  auto &irrep = block.irrep();

  // get limits of up
  auto it_up = block.ups_lower_upper_.find(up);
  auto [dn_lower, dn_upper] = it_up->second;

  // Get limits of flipped and represeantative up
  bit_t up_flip = up ^ mask;
  auto [up_rep, n_stable_syms, stable_syms] =
      group_action.representative_indices(up_flip);
  auto it_flip = block.ups_lower_upper_.find(up_rep);

  if (it_flip == block.ups_lower_upper_.end()) {
    return;
  }
  auto [dn_lower_rep, dn_upper_rep] = it_flip->second;

  // Loop over all dn configurations for this up configuration
  for (idx_t idx_in = dn_lower; idx_in < dn_upper; ++idx_in) {
    bit_t dn = block.dns_[idx_in];

    if (dn & dnmask) {
      // Compute flipped representative and symmetry leading to it
      bit_t dn_flip = dn ^ mask;
      bit_t dn_rep = group_action.apply(stable_syms[0], dn_flip);
      int rep_sym = stable_syms[0];

      // loop over stable symmetries if stabilizer exists on ups
      if (n_stable_syms > 1) {
        for (int n_sym = 1; n_sym < n_stable_syms; ++n_sym) {
          bit_t dn_trans = group_action.apply(stable_syms[n_sym], dn_flip);
          if (dn_trans < dn_rep) {
            dn_rep = dn_trans;
            rep_sym = stable_syms[n_sym];
          }
        }
      }

      // find the index of the representative
      auto it = std::lower_bound(block.dns_.begin() + dn_lower_rep,
                                 block.dns_.begin() + dn_upper_rep, dn_rep);
      // dn_rep not found, this is due to it having zero norm
      if ((it == block.dns_.begin() + dn_upper_rep) || (*it != dn_rep)) {
        continue;
      }
      idx_t idx_out = std::distance(block.dns_.begin(), it);

      // Compute matrix element
      double fermi_up = group_action.fermi_sign(rep_sym, up_flip);
      double fermi_dn = group_action.fermi_sign(rep_sym, dn_flip);
      coeff_t val = Jhalf * fermi_up * fermi_dn *
                    complex_to<coeff_t>(irrep.character(rep_sym)) *
                    block.norm(idx_out) / block.norm(idx_in);

      // Fermi sign for dn
      if (popcnt(dn & spacemask) & 1)
        fill(idx_out, idx_in, -val);
      else
        fill(idx_out, idx_in, val);
    }
  }
}

template <class bit_t, class coeff_t, class GroupAction, class Filler>
void do_exchange_symmetric(BondList const &bonds, Couplings const &couplings,
                           tJSymmetricSimple<bit_t, GroupAction> const &block,
                           Filler &&fill) {
  using bitops::gbit;
  using bitops::popcnt;

  auto exchange = bonds.bonds_of_type("HEISENBERG") +
                  bonds.bonds_of_type("EXCHANGE") + bonds.bonds_of_type("HB");
  auto exchange_tj =
      bonds.bonds_of_type("TJHEISENBERG") + bonds.bonds_of_type("TJHB");

  for (auto bond : exchange + exchange_tj) {

    if (bond.size() != 2)
      lila::Log.err("Error computing tJ exchange: "
                    "bonds must have exactly two sites defined");

    if (!utils::coupling_is_zero(bond, couplings)) {
      int s1 = std::min(bond.site(0), bond.site(1));
      int s2 = std::max(bond.site(0), bond.site(1));
      bit_t s1mask = (bit_t)1 << s1;
      bit_t s2mask = (bit_t)1 << s2;
      bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
      bit_t spacemask = ((1 << (s2 - s1 - 1)) - 1)
                        << (s1 + 1); // used to determine Fermi sign

      std::string cpl = bond.coupling();
      double J = lila::real(couplings[cpl]);
      double Jhalf = J / 2.;

      idx_t idx_up = 0;

      // Loop over all up configurations
      for (auto [up, lower_upper] : block.ups_lower_upper_) {
        if ((popcnt(up & mask) == 2) || (popcnt(up & mask) == 0))
          continue;

        // lower s1, raise, s2
        if ((up & mask) == s1mask) {

          // decide Fermi sign of upspins
          if (popcnt(up & spacemask) & 1)
            do_down_flips<bit_t, coeff_t>(up, idx_up, mask, spacemask, s2mask,
                                          Jhalf, block, fill);
          else
            do_down_flips<bit_t, coeff_t>(up, idx_up, mask, spacemask, s2mask,
                                          -Jhalf, block, fill);
          // lower s2, raise, s1
        } else if ((up & mask) == s2mask) {

          // decide Fermi sign of upspins
          if (popcnt(up & spacemask) & 1)
            do_down_flips<bit_t, coeff_t>(up, idx_up, mask, spacemask, s1mask,
                                          Jhalf, block, fill);
          else
            do_down_flips<bit_t, coeff_t>(up, idx_up, mask, spacemask, s1mask,
                                          -Jhalf, block, fill);
        }
        ++idx_up;
      }
    }
  }
}
} // namespace hydra::terms::tj_symmetric_simple
