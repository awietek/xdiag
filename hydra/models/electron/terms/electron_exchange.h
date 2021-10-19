#pragma once

#include <lila/utils/logger.h>

#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/models/electron/electron.h>
#include <hydra/models/utils/model_utils.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>

namespace hydra::electronterms {

template <class bit_t, class Filler>
void do_down_flips(bit_t up, idx_t idx_up, bit_t flipmask, bit_t spacemask,
                   bit_t dnmask, double val, Electron<bit_t> const &block,
                   Filler &&fill) {
  using utils::popcnt;

  int n_sites = block.n_sites();
  int n_dn = block.n_dn();

  // Get limits of flipped up
  bit_t up_flip = up ^ flipmask;
  idx_t idx_up_flip = block.lintable_up().index(up_flip);

  idx_t idx_out_offset = idx_up_flip * combinatorics::binomial(n_sites, n_dn);
  idx_t idx_in = idx_up * combinatorics::binomial(n_sites, n_dn);

  for (auto dn : Combinations(n_sites, n_dn)) {

    if ((dn & flipmask) == dnmask) {
      bit_t dn_flip = dn ^ flipmask;
      idx_t idx_out = idx_out_offset + block.lintable_dn().index(dn_flip);

      // decide Fermi sign from down spins
      if (popcnt(dn & spacemask) & 1)
        fill(idx_out, idx_in, -val);
      else
        fill(idx_out, idx_in, val);
    }

    ++idx_in;
  }
}

template <class bit_t, class Filler>
void do_exchange(BondList const &bonds, Couplings const &couplings,
                 Electron<bit_t> const &block, Filler &&fill) {
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
      lila::Log.err("Error computing Electron Exchange: "
                    "bond must have exactly two sites defined");

    std::string cpl = bond.coupling();
    if (!utils::coupling_is_zero(bond, couplings)) {
      int s1 = std::min(bond.site(0), bond.site(1));
      int s2 = std::max(bond.site(0), bond.site(1));
      bit_t s1mask = (bit_t)1 << s1;
      bit_t s2mask = (bit_t)1 << s2;
      bit_t flipmask = s1mask | s2mask;
      bit_t spacemask = ((1 << (s2 - s1 - 1)) - 1) << (s1 + 1);

      double J = lila::real(couplings[cpl]);
      double val = J / 2.;

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

} // namespace hydra::electronterms
