#pragma once
#include <hydra/common.h>

#include <hydra/blocks/spinhalf_symmetric/spinhalf_symmetric.h>
#include <hydra/blocks/utils/block_utils.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <hydra/bitops/bitops.h>

namespace hydra::terms::spinhalf_symmetric {

template <typename bit_t, typename coeff_t, class Filler>
void do_ising_symmetric(
    BondList const &bonds, Couplings const &couplings,
    indexing::SpinhalfSymmetricIndexing<bit_t> const &indexing, Filler &&fill) {
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

    for (idx_t idx = 0; idx < indexing.size(); ++idx) {
      bit_t state = indexing.state(idx);
      if (bitops::popcnt(state & mask) & 1) {
        fill(idx, idx, val_diff);
      } else {
        fill(idx, idx, val_same);
      }
    }
  }
}

} // namespace hydra::terms::spinhalf_symmetric
