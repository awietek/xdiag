#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/model_utils.h>
#include <hydra/symmetries/symmetry_utils.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::tj {

template <class bit_t, class GroupAction, class Filler>
void do_ising_symmetric(BondList const &bonds, Couplings const &couplings,
                        tJSymmetric<bit_t, GroupAction> const &block,
                        Filler &&fill) {
  using utils::gbit;
  using utils::popcnt;

  auto ising = bonds.bonds_of_type("HEISENBERG") +
               bonds.bonds_of_type("ISING") + bonds.bonds_of_type("HB");
  auto ising_tj = bonds.bonds_of_type("TJHEISENBERG") +
                  bonds.bonds_of_type("TJISING") + bonds.bonds_of_type("TJHB");

  for (auto bond : ising + ising_tj) {

    if (bond.size() != 2)
      lila::Log.err("Error computing tJ Ising: "
                    "bond must have exactly two sites defined");

    if (!utils::coupling_is_zero(bond, couplings)) {

      double J = lila::real(couplings[bond.coupling()]);

      // Set values for same/diff (tJ model definition)
      std::string type = bond.type();
      double val_same, val_diff;
      if ((type == "HEISENBERG") || (type == "ISING") || (type == "HB")) {
        val_same = J / 4.;
        val_diff = -J / 4.;
      } else if ((type == "TJHEISENBERG") || (type == "TJISING") ||
                 (type == "TJHB")) {
        val_same = 0.;
        val_diff = -J / 2.;
      }

      int s1 = bond.site(0);
      int s2 = bond.site(1);
      bit_t s1mask = (bit_t)1 << s1;
      bit_t s2mask = (bit_t)1 << s2;
      bit_t mask = s1mask | s2mask;

      for (auto [up, lower_upper] : block.ups_lower_upper_) {
        idx_t lower = lower_upper.first;
        idx_t upper = lower_upper.second;

        if (popcnt(up & mask) == 2) { // both spins pointing up
          for (idx_t idx = lower; idx < upper; ++idx) {
            fill(idx, idx, val_same);
          }
        } else if (up & s1mask) { // s1 is pointing up
          for (idx_t idx = lower; idx < upper; ++idx) {
            bit_t dn = block.dns_[idx];
            if (dn & s2mask)
              fill(idx, idx, val_diff);
          }
        } else if (up & s2mask) { // s2 is pointing up
          for (idx_t idx = lower; idx < upper; ++idx) {
            bit_t dn = block.dns_[idx];
            if (dn & s1mask)
              fill(idx, idx, val_diff);
          }
        } else { // no upspins
          for (idx_t idx = lower; idx < upper; ++idx) {
            bit_t dn = block.dns_[idx];
            if (popcnt(dn & mask) == 2)
              fill(idx, idx, val_same);
          }
        }
      }
    }
  }
}
} // namespace hydra::tj
