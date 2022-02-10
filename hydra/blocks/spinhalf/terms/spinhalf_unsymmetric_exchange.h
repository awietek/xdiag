#pragma once

#include <hydra/common.h>

#include <hydra/bitops/bitops.h>

#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/utils/block_utils.h>

#include <hydra/combinatorics/combinations.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/operators/operator_utils.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void spinhalf_unsymmetric_exchange(BondList const &bonds,
                                   Couplings const &couplings,
                                   Indexing &&indexing, Filler &&fill) {
  using bitops::gbit;

  auto clean_bonds = utils::clean_bondlist(bonds, couplings,
                                           {"HEISENBERG", "HB", "EXCHANGE"}, 2);
  for (auto bond : clean_bonds) {
    std::string cpl = bond.coupling();

    utils::check_sites_disjoint(bond);
    int s1 = bond[0];
    int s2 = bond[1];
    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

    auto [J, J_conj] = utils::get_coupling_and_conj<coeff_t>(couplings, cpl);
    coeff_t Jhalf = J / 2.;
    coeff_t Jhalf_conj = J_conj / 2.;

    idx_t idx_in = 0;
    for (bit_t spins : indexing.states()) {
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

} // namespace hydra::terms
