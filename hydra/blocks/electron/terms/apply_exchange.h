#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>
#include <hydra/utils/print_macro.h>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void electron_do_down_flips(bit_t ups, idx_t idx_ups, bit_t flipmask,
                            bit_t fermimask, int sdn, coeff_t val,
                            Indexing &&indexing, Filler &&fill) {
  idx_t size_dns = indexing.size_dns();

  // Get limits of flipped up
  bit_t ups_flip = ups ^ flipmask;
  idx_t idx_ups_flip = indexing.index_ups(ups_flip);

  idx_t idx_out_offset = idx_ups_flip * size_dns;
  idx_t idx_in = idx_ups * size_dns;

  for (auto dns : indexing.states_dns()) {

    if ((bitops::popcnt(dns & flipmask) == 1) &&
        ((bool)bitops::gbit(dns, sdn))) {
      bit_t dns_flip = dns ^ flipmask;
      idx_t idx_out = idx_out_offset + indexing.index_dns(dns_flip);
      fill(idx_out, idx_in, (bitops::popcnt(dns & fermimask) & 1) ? -val : val);
    }

    ++idx_in;
  }
}

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void apply_exchange(Bond const &bond, Indexing &&indexing, Filler &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined() && (bond.type() == "EXCHANGE"));
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  coeff_t J = bond.coupling<coeff_t>();
  int s1 = std::min(bond.site(0), bond.site(1));
  int s2 = std::max(bond.site(0), bond.site(1));

  // Prepare bitmasks
  bit_t s1mask = (bit_t)1 << s1;
  bit_t s2mask = (bit_t)1 << s2;
  bit_t flipmask = s1mask | s2mask;
  bit_t fermimask = ((1 << (s2 - s1 - 1)) - 1) << (s1 + 1);

  idx_t idx_up = 0;
  for (auto ups : indexing.states_ups()) {

    if (bitops::popcnt(ups & flipmask) == 1) {

      // Set the correct prefactor
      coeff_t Jhalf;
      if constexpr (is_complex<coeff_t>()) {
        Jhalf = (bitops::gbit(ups, s1)) ? J / 2.0 : hydra::conj(J) / 2.0;
      } else {
        Jhalf = J / 2.;
      }

      // decide Fermi sign of upspins
      if (bitops::popcnt(ups & fermimask) & 1) {
        electron_do_down_flips(ups, idx_up, flipmask, fermimask,
                               bitops::gbit(ups, s1) ? s2 : s1, Jhalf, indexing,
                               fill);
      } else {
        electron_do_down_flips(ups, idx_up, flipmask, fermimask,
                               bitops::gbit(ups, s1) ? s2 : s1, -Jhalf,
                               indexing, fill);
      }
    }
    ++idx_up;
  }
}

} // namespace hydra::electron
