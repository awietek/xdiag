#pragma once

#include <hydra/common.h>

#include <hydra/bitops/bitops.h>

#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/utils/block_utils.h>

#include <hydra/combinatorics/combinations.h>

#include <hydra/indexing/tj/tj_indexing.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/operators/operator_utils.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void tj_ising(BondList const &bonds, Couplings const &couplings,
              Indexing &&indexing, Filler &&fill) {
  using bitops::gbit;
  using bitops::popcnt;
  using combinatorics::Combinations;

  int n_sites = indexing.n_sites();
  int n_up = indexing.n_up();
  int n_holes = indexing.n_holes();
  int n_spins = indexing.n_spins();
  idx_t size_spins = indexing.size_spins();

  auto clean_bonds = utils::clean_bondlist(
      bonds, couplings,
      {"HEISENBERG", "HB", "ISING", "TJHEISENBERG", "TJISING", "TJHB"}, 2);

  for (auto bond : clean_bonds) {

    std::string type = bond.type();
    std::string cpl = bond.coupling();

    utils::check_sites_disjoint(bond);
    int s1 = bond[0];
    int s2 = bond[1];

    coeff_t J = utils::get_coupling<coeff_t>(couplings, cpl);

    // Set values for same/diff (tJ block definition)
    coeff_t val_same, val_diff;
    if ((type == "HEISENBERG") || (type == "ISING") || (type == "HB")) {
      val_same = J / 4.;
      val_diff = -J / 4.;
    } else if ((type == "TJHEISENBERG") || (type == "TJISING") ||
               (type == "TJHB")) {
      val_same = 0.;
      val_diff = -J / 2.;
    } else {
      continue;
    }

    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    bit_t sitesmask = ((bit_t)1 << n_sites) - 1;
    // bit_t spinsmask = ((bit_t)1 << n_spins) - 1;

    idx_t idx = 0;
    for (auto holes : Combinations<bit_t>(n_sites, n_holes)) {

      // there cannot be a hole on one site
      if (popcnt(holes & flipmask) == 0) {

        bit_t not_holes = (~holes) & sitesmask;
        for (auto spins : Combinations<bit_t>(n_spins, n_up)) {
          bit_t ups_dns = bitops::deposit(spins, not_holes);

          if (gbit(ups_dns, s1) == gbit(ups_dns, s2)) {
            fill(idx, idx, val_same);
          } else {
            fill(idx, idx, val_diff);
          }
          ++idx;
        }
        // a hole is present at site -> no Ising term -> skip ahead
      } else {
        idx += size_spins;
      }
    }
  }
}

} // namespace hydra::terms
