#pragma once

#include <hydra/common.h>

#include <hydra/bitops/bitops.h>

#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/utils/block_utils.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void spinhalf_symmetric_exchange(BondList const &bonds,
                                 Couplings const &couplings,
                                 Indexing &&indexing, Filler &&fill) {
  using bitops::gbit;

  auto bloch_factors = utils::characters<coeff_t>(indexing.irrep());
  auto clean_bonds = utils::clean_bondlist(bonds, couplings,
                                           {"HEISENBERG", "HB", "EXCHANGE"}, 2);
  for (auto bond : clean_bonds) {

    // Get sites and define flipmask
    utils::check_sites_disjoint(bond);
    int s1 = bond[0];
    int s2 = bond[1];
    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

    // Get couplings
    std::string cpl = bond.coupling();
    auto [J, J_conj] = utils::get_coupling_and_conj<coeff_t>(couplings, cpl);
    coeff_t Jhalf = J / 2.;
    coeff_t Jhalf_conj = J_conj / 2.;

    // Run through all states and apply bond
    for (idx_t idx_in = 0; idx_in < indexing.size(); ++idx_in) {
      bit_t spins_in = indexing.state(idx_in);

      if (bitops::popcnt(spins_in & flipmask) & 1) { // if spins are flippable

        bit_t spins_out = spins_in ^ flipmask;
        auto [idx_out, syms] = indexing.index_syms(spins_out);

        // if new spins has non-zero norm, compute element and fill
        if (idx_out != invalid_index) {
          double norm_out = indexing.norm(idx_out);
          double norm_in = indexing.norm(idx_in);
          coeff_t bloch = bloch_factors[syms[0]];
          if constexpr (is_complex<coeff_t>()) {
            coeff_t val = (gbit(spins_in, s1) ? Jhalf : Jhalf_conj) * bloch *
                          norm_out / norm_in;
            fill(idx_out, idx_in, val);
          } else {
            coeff_t val = Jhalf * bloch * norm_out / norm_in;
            fill(idx_out, idx_in, val);
          }
        }
      }
    }
  }
}

} // namespace hydra::terms
