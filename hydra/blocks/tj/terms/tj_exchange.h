#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/tj/tj_utils.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/operators/operator_utils.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void tj_exchange(BondList const &bonds, Couplings const &couplings,
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
      {"HEISENBERG", "HB", "EXCHANGE", "TJHEISENBERG", "TJISING", "TJHB"}, 2);

  for (auto bond : clean_bonds) {

    std::string type = bond.type();
    std::string cpl = bond.coupling();

    auto [J, J_conj] = utils::get_coupling_and_conj<coeff_t>(couplings, cpl);
    coeff_t Jhalf = J / 2.;
    coeff_t Jhalf_conj = J_conj / 2.;
    if constexpr (is_real<coeff_t>()) {  // just to supress "unused" warning
       assert(Jhalf == Jhalf_conj);
    }
    
    utils::check_sites_disjoint(bond);
    int s1 = bond[0];
    int s2 = bond[1];

    // Prepare bitmasks
    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    bit_t sitesmask = ((bit_t)1 << n_sites) - 1;
    int l = std::min(s1, s2);
    int u = std::max(s1, s2);
    bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

    // bit_t spinsmask = ((bit_t)1 << n_spins) - 1;

    idx_t idx = 0;
    idx_t holes_idx = 0;
    for (auto holes : Combinations<bit_t>(n_sites, n_holes)) {

      // no holes should be present for exchange
      if (popcnt(holes & flipmask) == 0) {

        bit_t not_holes = (~holes) & sitesmask;
        idx_t holes_offset = holes_idx * size_spins;

        for (auto spins : Combinations<bit_t>(n_spins, n_up)) {
          bit_t ups = bitops::deposit(spins, not_holes);

          if (popcnt(ups & flipmask) & 1) { // spins are flippable
            bit_t new_ups = ups ^ flipmask;
            bit_t new_spins = bitops::extract(new_ups, not_holes);
            idx_t new_idx = holes_offset + indexing.index_spins(new_spins);

            bool fermi_sign = popcnt(not_holes & fermimask) & 1;

            if constexpr (is_complex<coeff_t>()) {
              if (gbit(ups, s1)) {
                fill(new_idx, idx, (fermi_sign) ? Jhalf : -Jhalf);
              } else {
                fill(new_idx, idx, (fermi_sign) ? Jhalf_conj : -Jhalf_conj);
              }
            } else {
              fill(new_idx, idx, (fermi_sign) ? Jhalf : -Jhalf);
            }
          }
          ++idx;
        }
      } else {
        idx += size_spins;
      }
      ++holes_idx;
    } // for (auto holes : ...)

  } // for (auto bond : clean_bonds)
}

} // namespace hydra::terms
