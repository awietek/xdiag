#pragma once
#include <hydra/common.h>

#include <hydra/combinatorics/combinations.h>

#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/utils/block_utils.h>

#include <hydra/indexing/spinhalf/spinhalf_indexing.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <hydra/bitops/bitops.h>

namespace hydra::terms::spinhalf {

template <class bit_t, class coeff_t, class Filler>
void do_ising(BondList const &bonds, Couplings const &couplings,
              indexing::SpinhalfIndexing<bit_t> const &indexing,
              Filler &&fill) {
  using combinatorics::Combinations;

  int n_sites = indexing.n_sites();
  int n_up = indexing.n_up();

  auto clean_bonds =
      utils::clean_bondlist(bonds, couplings, {"HEISENBERG", "HB", "ISING"}, 2);
  for (auto bond : clean_bonds) {

    utils::check_sites_disjoint(bond);
    int s1 = bond[0];
    int s2 = bond[1];
    bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

    std::string cpl = bond.coupling();
    coeff_t J = utils::get_coupling<coeff_t>(couplings, cpl);
    coeff_t val_same = J / 4.;
    coeff_t val_diff = -J / 4.;

    idx_t idx = 0;
    for (auto spins : Combinations<bit_t>(n_sites, n_up)) {

      if (bitops::popcnt(spins & mask) & 1)
        fill(idx, idx, val_diff);
      else
        fill(idx, idx, val_same);

      ++idx;
    }
  }
}

} // namespace hydra::terms::spinhalf
