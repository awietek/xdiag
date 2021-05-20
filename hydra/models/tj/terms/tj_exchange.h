#pragma once

#include <hydra/common.h>
#include <hydra/models/tj/tj.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>

namespace hydra::tjdetail {

template <class bit_t, class Filler>
void do_down_flips(bit_t up, idx_t idx_up, bit_t up_mask, bit_t dn_mask,
                   double val, tJ<bit_t> const &block, Filler &&fill) {
  bit_t mask = up_mask | dn_mask;



  // Get limits of flipped up
  bit_t up_flip = up ^ mask;
  idx_t idx_up_flip = block.lintable_up_.index(up_flip);
  auto [dn_lower_flip, dn_upper_flip] = block.dn_limits_for_up_[idx_up_flip];

  // get limits of up
  auto [dn_lower, dn_upper] = block.dn_limits_for_up_[idx_up];

  for (idx_t idx_in = dn_lower; idx_in < dn_upper; ++idx_in) {
    bit_t dn = block.dns_[idx_in];

    if (dn & dn_mask) {
      bit_t dn_flip = dn ^ mask;
      auto it = std::lower_bound(block.dns_.begin() + dn_lower_flip,
                                 block.dns_.begin() + dn_upper_flip, dn_flip);
      idx_t idx_out = std::distance(block.dns_.begin(), it);
      fill(idx_out, idx_in, val);
    }
  }
}

template <class bit_t, class Filler>
void do_exchange(BondList const &bonds, Couplings const &couplings,
                 tJ<bit_t> const &block, Filler &&fill) {
  using utils::gbit;
  using utils::popcnt;

  int n_sites = block.n_sites();
  int n_up = block.n_up();

  auto exchange = bonds.bonds_of_type("HEISENBERG") +
                  bonds.bonds_of_type("EXCHANGE") + bonds.bonds_of_type("HB");
  auto exchange_tj =
      bonds.bonds_of_type("TJHEISENBERG") + bonds.bonds_of_type("TJHB");

  for (auto bond : exchange + exchange_tj) {

    if (bond.size() != 2)
      HydraLog.err("Error computing tJ Exchange: "
                   "bond must have exactly two sites defined");

    std::string cpl = bond.coupling();
    if (couplings.defined(cpl) && !lila::close(couplings[cpl], (complex)0.)) {

      double J = lila::real(couplings[cpl]);
      double val = J / 2.;

      int s1 = bond.site(0);
      int s2 = bond.site(1);
      bit_t s1mask = (bit_t)1 << s1;
      bit_t s2mask = (bit_t)1 << s2;
      bit_t mask = s1mask | s2mask;

      idx_t idx_up = 0;
      for (auto up : Combinations<bit_t>(n_sites, n_up)) {

        if ((up & mask) == s1mask)
          do_down_flips(up, idx_up, s1mask, s2mask, val, block, fill);
        else if ((up & mask) == s2mask)
          do_down_flips(up, idx_up, s2mask, s1mask, val, block, fill);

        ++idx_up;
      }
    }
  }
}

} // namespace hydra::tjdetail
