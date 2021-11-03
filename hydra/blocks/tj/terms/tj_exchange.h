#pragma once

#include <lila/utils/logger.h>

#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/bitops/bitops.h>

namespace hydra::terms::tj {

template <class bit_t, class Filler>
void do_down_flips(bit_t up, idx_t idx_up, bit_t flipmask, bit_t spacemask,
                   bit_t dnmask, double val, tJ<bit_t> const &block,
                   Filler &&fill) {
  using bitops::popcnt;

  // Get limits of flipped up
  bit_t up_flip = up ^ flipmask;
  idx_t idx_up_flip = block.lintable_up_.index(up_flip);
  auto [dn_lower_flip, dn_upper_flip] = block.dn_limits_for_up_[idx_up_flip];

  // get limits of up
  auto [dn_lower, dn_upper] = block.dn_limits_for_up_[idx_up];

  for (idx_t idx_in = dn_lower; idx_in < dn_upper; ++idx_in) {
    bit_t dn = block.dns_[idx_in];

    if (dn & dnmask) {
      bit_t dn_flip = dn ^ flipmask;
      auto it = std::lower_bound(block.dns_.begin() + dn_lower_flip,
                                 block.dns_.begin() + dn_upper_flip, dn_flip);
      idx_t idx_out = std::distance(block.dns_.begin(), it);

      // decide Fermi sign from down spins
      if (popcnt(dn & spacemask) & 1)
        fill(idx_out, idx_in, -val);
      else
        fill(idx_out, idx_in, val);
    }
  }
}

template <class bit_t, class Filler>
void do_exchange(BondList const &bonds, Couplings const &couplings,
                 tJ<bit_t> const &block, Filler &&fill) {
  using bitops::gbit;
  using bitops::popcnt;
  using combinatorics::Combinations;

  int n_sites = block.n_sites();
  int n_up = block.n_up();

  auto exchange = bonds.bonds_of_type("HEISENBERG") +
                  bonds.bonds_of_type("EXCHANGE") + bonds.bonds_of_type("HB");
  auto exchange_tj =
      bonds.bonds_of_type("TJHEISENBERG") + bonds.bonds_of_type("TJHB");

  for (auto bond : exchange + exchange_tj) {

    if (bond.size() != 2)
      lila::Log.err("Error computing tJ Exchange: "
                    "bond must have exactly two sites defined");

    std::string cpl = bond.coupling();
    if (!utils::coupling_is_zero(bond, couplings)) {

      double J = lila::real(couplings[cpl]);
      double val = J / 2.;

      int s1 = std::min(bond.site(0), bond.site(1));
      int s2 = std::max(bond.site(0), bond.site(1));
      bit_t s1mask = (bit_t)1 << s1;
      bit_t s2mask = (bit_t)1 << s2;
      bit_t flipmask = s1mask | s2mask;
      bit_t spacemask = ((1 << (s2 - s1 - 1)) - 1) << (s1 + 1);

      idx_t idx_up = 0;
      for (auto up : Combinations<bit_t>(n_sites, n_up)) {

        // lower s1, raise, s2
        if ((up & flipmask) == s1mask) {

          // decide Fermi sign of upspins
          if (popcnt(up & spacemask) & 1)
            do_down_flips(up, idx_up, flipmask, spacemask, s2mask, val, block,
                          fill);
          else
            do_down_flips(up, idx_up, flipmask, spacemask, s2mask, -val, block,
                          fill);
          // lower s2, raise, s1
        } else if ((up & flipmask) == s2mask) {

          // decide Fermi sign of upspins
          if (popcnt(up & spacemask) & 1)
            do_down_flips(up, idx_up, flipmask, spacemask, s1mask, val, block,
                          fill);
          else
            do_down_flips(up, idx_up, flipmask, spacemask, s1mask, -val, block,
                          fill);
        }

        ++idx_up;
      }
    }
  }
}

} // namespace hydra::tjdetail
