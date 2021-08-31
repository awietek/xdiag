#pragma once

#include <lila/utils/logger.h>

#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/models/tj/tj.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>

namespace hydra::tjterms {

template <class bit_t, class Filler>
void do_ising(BondList const &bonds, Couplings const &couplings,
              tJ<bit_t> const &block, Filler &&fill) {
  using utils::gbit;
  using utils::popcnt;

  int n_sites = block.n_sites();
  int n_up = block.n_up();

  auto ising = bonds.bonds_of_type("HEISENBERG") +
               bonds.bonds_of_type("ISING") + bonds.bonds_of_type("HB");
  auto ising_tj = bonds.bonds_of_type("TJHEISENBERG") +
                  bonds.bonds_of_type("TJISING") + bonds.bonds_of_type("TJHB");

  for (auto bond : ising + ising_tj) {

    if (bond.size() != 2)
      lila::Log.err("Error computing tJ Ising: "
                    "bond must have exactly two sites defined");

    std::string coupling = bond.coupling();
    if (couplings.defined(coupling) &&
        !lila::close(couplings[coupling], (complex)0.)) {

      double J = lila::real(couplings[coupling]);

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

      idx_t idx_up = 0;
      for (auto up : Combinations<bit_t>(n_sites, n_up)) {
        auto [dn_lower, dn_upper] = block.dn_limits_for_up_[idx_up];

        if (popcnt(up & mask) == 2) { // both spins pointing up
          for (idx_t idx = dn_lower; idx < dn_upper; ++idx) {
            fill(idx, idx, val_same);
          }
        } else if (up & s1mask) { // s1 is pointing up
          for (idx_t idx = dn_lower; idx < dn_upper; ++idx) {
            bit_t dn = block.dns_[idx];
            if (dn & s2mask)
              fill(idx, idx, val_diff);
          }
        } else if (up & s2mask) { // s2 is pointing up
          for (idx_t idx = dn_lower; idx < dn_upper; ++idx) {
            bit_t dn = block.dns_[idx];
            if (dn & s1mask)
              fill(idx, idx, val_diff);
          }
        } else { // no upspins
          for (idx_t idx = dn_lower; idx < dn_upper; ++idx) {
            bit_t dn = block.dns_[idx];
            if (popcnt(dn & mask) == 2)
              fill(idx, idx, val_same);
          }
        }

        ++idx_up;
      }
    }
  }
}

} // namespace hydra::tjdetail
