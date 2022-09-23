#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void apply_ising(Bond const &bond, Indexing &&indexing, Filler &&fill) {
  using combinatorics::Combinations;

  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert((bond.type() == "ISING") || (bond.type() == "TJISING"));
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  int n_sites = indexing.n_sites();
  int n_up = indexing.n_up();
  int n_holes = indexing.n_holes();
  int n_spins = indexing.n_spins();
  idx_t size_spins = indexing.size_spins();

  int s1 = bond[0];
  int s2 = bond[1];
  coeff_t J = bond.coupling<coeff_t>();

  // Set values for same/diff (tJ block definition)
  coeff_t val_same, val_diff;
  if (bond.type() == "ISING") {
    val_same = J / 4.;
    val_diff = -J / 4.;
  } else {  // (bond.type() == "TJISING")
    val_same = 0.;
    val_diff = -J / 2.;
  }

  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;
  // bit_t spinsmask = ((bit_t)1 << n_spins) - 1;

  idx_t idx = 0;
  for (auto holes : Combinations<bit_t>(n_sites, n_holes)) {

    // there cannot be a hole on one site
    if (bitops::popcnt(holes & flipmask) == 0) {

      bit_t not_holes = (~holes) & sitesmask;
      for (auto spins : Combinations<bit_t>(n_spins, n_up)) {
        bit_t ups_dns = bitops::deposit(spins, not_holes);

        if (bitops::gbit(ups_dns, s1) == bitops::gbit(ups_dns, s2)) {
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

} // namespace hydra::tj
