#pragma once

#include <lila/utils/logger.h>

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/indexing/spinhalf/spinhalf_symmetric_indexing.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra::terms::spinhalf {

template <typename bit_t, typename coeff_t, typename Filler>
void do_exchange(BondList const &bonds, Couplings const &couplings,
                 indexing::SpinhalfIndexing<bit_t> const &indexing,
                 Filler &&fill) {
  using bitops::gbit;
  using combinatorics::Combinations;

  int n_sites = indexing.n_sites();
  int n_up = indexing.n_up();

  auto clean_bonds = utils::clean_bondlist(bonds, couplings,
                                           {"HEISENBERG", "HB", "EXCHANGE"}, 2);
  for (auto bond : clean_bonds) {

    int s1 = bond[0];
    int s2 = bond[1];
    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

    if (s1 == s2)
      lila::Log.err("Error: Exchange operator acting on twice the same site");

    std::string cpl = bond.coupling();
    auto [J, J_conj] = utils::get_coupling_and_conj<coeff_t>(couplings, cpl);
    coeff_t Jhalf = J / 2.;
    coeff_t Jhalf_conj = J_conj / 2.;

    idx_t idx_in = 0;
    for (bit_t spins : Combinations<bit_t>(n_sites, n_up)) {
      if (bitops::popcnt(spins & flipmask) & 1) {
        bit_t new_spins = spins ^ flipmask;
        idx_t idx_out = indexing.index(new_spins);

        if constexpr (is_complex<coeff_t>()) {
          fill(idx_out, idx_in, gbit(spins, s1) ? Jhalf : Jhalf_conj);
        } else {
          fill(idx_out, idx_in, Jhalf);
        }
      }
      ++idx_in;
    }
  }
}

} // namespace hydra::terms::spinhalf
