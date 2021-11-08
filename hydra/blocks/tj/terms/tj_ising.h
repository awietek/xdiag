#pragma once

#include <lila/utils/logger.h>

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/tj/tj_utils.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms::tj {

template <class bit_t, class Filler>
void do_ising(BondList const &bonds, Couplings const &couplings,
              tJ<bit_t> const &block, Filler &&fill) {
  using bitops::bits_to_string;
  using bitops::gbit;
  using bitops::popcnt;
  using combinatorics::Combinations;

  int n_sites = block.n_sites();
  int nup = block.n_up();
  int ndn = block.n_dn();
  int charge = nup + ndn;
  int n_holes = n_sites - charge;
  idx_t size_spins = block.size_spins_;

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

      int s1 = bond[0];
      int s2 = bond[1];
      bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
      bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

      idx_t idx = 0;
      for (auto holes : Combinations<bit_t>(n_sites, n_holes)) {

        // a hole is present at site -> no Ising term
        if (popcnt(holes & flipmask) != 0) {
          idx += size_spins;
          continue;
        }

        bit_t not_holes = (~holes) & sitesmask;
        for (auto spins : Combinations<bit_t>(charge, nup)) {
          bit_t ups_dns = bitops::deposit(spins, not_holes);
          if (gbit(ups_dns, s1) == gbit(ups_dns, s2)) {
            fill(idx, idx, val_same);
          } else {
            fill(idx, idx, val_diff);
          }
          ++idx;
        }
      }
    }
  }
}

} // namespace hydra::terms::tj
