#pragma once

#include <lila/utils/logger.h>

#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/models/electron/electron.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>

namespace hydra::electronterms {

template <class bit_t, class Filler>
void do_ising(BondList const &bonds, Couplings const &couplings,
              Electron<bit_t> const &block, Filler &&fill) {
  using utils::gbit;
  using utils::popcnt;

  int n_sites = block.n_sites();
  int n_up = block.n_up();
  int n_dn = block.n_dn();

  auto ising = bonds.bonds_of_type("HEISENBERG") +
               bonds.bonds_of_type("ISING") + bonds.bonds_of_type("HB");
  auto ising_tj = bonds.bonds_of_type("TJHEISENBERG") +
                  bonds.bonds_of_type("TJISING") + bonds.bonds_of_type("TJHB");

  for (auto bond : ising + ising_tj) {

    if (bond.size() != 2)
      lila::Log.err("Error computing Electron Ising: "
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

      idx_t idx = 0;
      for (auto up : Combinations<bit_t>(n_sites, n_up)) {

        if ((up & mask) == mask) { // both spins pointing up
          for (bit_t dn : Combinations<bit_t>(n_sites, n_dn)) {
            if (!(dn & mask)) {
              fill(idx, idx, val_same);
            }
            ++idx;
          }
        } else if (up & s1mask) { // s1 is pointing up
          for (bit_t dn : Combinations<bit_t>(n_sites, n_dn)) {
            if ((dn & mask) == s2mask) {
              fill(idx, idx, val_diff);
            }
            ++idx;
          }
        } else if (up & s2mask) { // s2 is pointing up
          for (bit_t dn : Combinations<bit_t>(n_sites, n_dn)) {
            if ((dn & mask) == s1mask) {
              fill(idx, idx, val_diff);
            }
            ++idx;
          }

        } else { // no upspins
          for (bit_t dn : Combinations<bit_t>(n_sites, n_dn)) {
            if ((dn & mask) == mask) {
              fill(idx, idx, val_same);
            }
            ++idx;
          }
        }

      } // for (auto up : ...)
    }
  }
}

} // namespace hydra::electron
