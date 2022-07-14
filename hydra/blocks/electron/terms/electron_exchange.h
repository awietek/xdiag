#pragma once

#include <hydra/common.h>

#include <hydra/blocks/utils/block_utils.h>

#include <hydra/bitops/bitops.h>

#include <hydra/indexing/electron/electron_indexing.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/operators/operator_utils.h>

#include <hydra/combinatorics/combinations.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void electron_do_down_flips(bit_t ups, idx_t idx_ups, bit_t flipmask,
                            bit_t fermimask, int sdn, coeff_t val,
                            Indexing &&indexing, Filler &&fill) {
  using bitops::gbit;
  using bitops::popcnt;

  idx_t size_dns = indexing.size_dns();

  // Get limits of flipped up
  bit_t ups_flip = ups ^ flipmask;
  idx_t idx_ups_flip = indexing.index_ups(ups_flip);

  idx_t idx_out_offset = idx_ups_flip * size_dns;
  idx_t idx_in = idx_ups * size_dns;

  for (auto dns : indexing.states_dns()) {

    if ((popcnt(dns & flipmask) == 1) && ((bool)gbit(dns, sdn))) {
      bit_t dns_flip = dns ^ flipmask;
      idx_t idx_out = idx_out_offset + indexing.index_dns(dns_flip);
      fill(idx_out, idx_in, (popcnt(dns & fermimask) & 1) ? -val : val);
    }

    ++idx_in;
  }
}

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void electron_exchange(BondList const &bonds, Couplings const &couplings,
                       Indexing &&indexing, Filler &&fill) {

  using bitops::gbit;
  using bitops::popcnt;

  auto clean_bonds = utils::clean_bondlist(
      bonds, couplings,
      {"HEISENBERG", "HB", "EXCHANGE", "TJHEISENBERG", "TJHB"}, 2);

  for (auto bond : clean_bonds) {
    std::string cpl = bond.coupling();

    utils::check_sites_disjoint(bond);
    int s1 = std::min(bond.site(0), bond.site(1));
    int s2 = std::max(bond.site(0), bond.site(1));

    auto [J, J_conj] = utils::get_coupling_and_conj<coeff_t>(couplings, cpl);
    if constexpr (is_real<coeff_t>()) {  // just to supress "unused" warning
	assert(J == J_conj);
      }
    
    // Prepare bitmasks
    bit_t s1mask = (bit_t)1 << s1;
    bit_t s2mask = (bit_t)1 << s2;
    bit_t flipmask = s1mask | s2mask;
    bit_t fermimask = ((1 << (s2 - s1 - 1)) - 1) << (s1 + 1);

    idx_t idx_up = 0;
    for (auto ups : indexing.states_ups()) {

      if (popcnt(ups & flipmask) == 1) {

        // Set the correct prefactor
        coeff_t Jhalf;
        if constexpr (is_complex<coeff_t>()) {
          Jhalf = (gbit(ups, s1)) ? J / 2. : J_conj / 2.;
        } else {
          Jhalf = J / 2.;
        }

        // decide Fermi sign of upspins
        if (popcnt(ups & fermimask) & 1) {
          electron_do_down_flips(ups, idx_up, flipmask, fermimask,
                                 gbit(ups, s1) ? s2 : s1, Jhalf, indexing,
                                 fill);
        } else {
          electron_do_down_flips(ups, idx_up, flipmask, fermimask,
                                 gbit(ups, s1) ? s2 : s1, -Jhalf, indexing,
                                 fill);
        }
      }
      ++idx_up;
    }
  }
}

} // namespace hydra::terms
