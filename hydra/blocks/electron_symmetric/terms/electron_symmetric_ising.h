#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/bitops/bitops.h>

#include <hydra/blocks/blocks.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/symmetry_utils.h>

namespace hydra::terms::electron_symmetric {

template <class bit_t, class GroupAction, class Filler>
void do_ising_symmetric(BondList const &bonds, Couplings const &couplings,
                        ElectronSymmetric<bit_t, GroupAction> const &block,
                        Filler &&fill) {

  using bitops::gbit;
  using bitops::popcnt;

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

      // Set values for same/diff (tJ block definition)
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

      idx_t idx = 0;
      for (bit_t ups : block.reps_up_) {
        auto const &dnss = block.dns_for_up_rep(ups);

        if ((ups & mask) == mask) { // both spins pointing up
          for (auto dns : dnss) {
            if (!(dns & mask))
              fill(idx, idx, val_same);
            ++idx;
          }
        } else if (ups & s1mask) { // s1 is pointing up
          for (auto dns : dnss) {
            if ((dns & mask) == s2mask)
              fill(idx, idx, val_diff);
            ++idx;
          }
        } else if (ups & s2mask) { // s2 is pointing up
          for (auto dns : dnss) {
            if ((dns & mask) == s1mask)
              fill(idx, idx, val_diff);
            ++idx;
          }
        } else { // no upspins
          for (auto dns : dnss) {
            if ((dns & mask) == mask)
              fill(idx, idx, val_same);
            ++idx;
          }
        }
      }
    }
  }
}

} // namespace hydra::terms::electron_symmetric
