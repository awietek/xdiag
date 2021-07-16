#pragma once

#include <lila/utils/logger.h>

#include <hydra/common.h>
#include <hydra/models/electron_symmetric/electron_symmetric.h>
#include <hydra/models/model_utils.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/symmetry_utils.h>
#include <hydra/utils/bitops.h>

namespace hydra::tj {

template <class bit_t, class coeff_t, class SymmetryGroup, class Filler>
void do_exchange_symmetric(BondList const &bonds, Couplings const &couplings,
                           tJSymmetric<bit_t, SymmetryGroup> const &block,
                           Filler &&fill) {
  using utils::gbit;
  using utils::popcnt;

  auto symmetry_group = block.symmetry_group();
  auto irrep = block.irrep();

  auto exchange = bonds.bonds_of_type("HEISENBERG") +
                  bonds.bonds_of_type("EXCHANGE") + bonds.bonds_of_type("HB");
  auto exchange_tj =
      bonds.bonds_of_type("TJHEISENBERG") + bonds.bonds_of_type("TJHB");

  for (auto bond : exchange + exchange_tj) {

    if (bond.size() != 2)
      lila::Log.err("Error computing Electron hopping: "
                    "hoppings must have exactly two sites defined");

    if (!utils::coupling_is_zero(bond, couplings)) {
      int s1 = bond.site(0);
      int s2 = bond.site(1);
      bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

      std::string cpl = bond.coupling();
      double J = lila::real(couplings[cpl]);
      double Jhalf = J / 2.;

      for (auto [up, lower_upper] : block.ups_lower_upper_) {
        if ((popcnt(up & mask) == 2) || (popcnt(up & mask) == 0))
          continue;

        bit_t up_flip = up ^ mask;
        auto [up_rep, n_stable_syms, stable_syms] =
            symmetry_group.representative_indices(up_flip);

        // Determine where to look for for dn configurations
        auto it1 = block.ups_lower_upper_.find(up_rep);
        auto [lower_flip, upper_flip] = it1->second;
        auto begin_flip = block.dns_.begin() + lower_flip;
        auto end_flip = block.dns_.begin() + upper_flip;

        // Loop over dn spins
        idx_t lower = lower_upper.first;
        idx_t upper = lower_upper.second;
        for (idx_t idx = lower; idx < upper; ++idx) {
          bit_t dn = block.dn(idx);
          if ((popcnt(dn & mask) != 1))
            continue;
          bit_t dn_flip = dn ^ mask;
          bit_t dn_rep = symmetry_group.apply(stable_syms[0], dn_flip);
          int rep_sym = stable_syms[0];

          // loop over stable symmetries if stabilizer exists on ups
          if (n_stable_syms > 1) {
            for (int n_sym = 1; n_sym < n_stable_syms; ++n_sym) {

              bit_t dn_trans =
                  symmetry_group.apply(stable_syms[n_sym], dn_flip);

              if (dn_trans < dn_rep) {
                dn_rep = dn_trans;
                rep_sym = stable_syms[n_sym];
              }
            }
          }

          auto it = std::lower_bound(begin_flip, end_flip, dn_rep);
          idx_t idx_out = std::distance(block.dns_.begin(), it);
          fill(idx_out, idx, Jhalf);
        }
      }
    }
  }
}
} // namespace hydra::tj

