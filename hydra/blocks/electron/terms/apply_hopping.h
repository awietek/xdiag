#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void apply_hopping(Bond const &bond, Indexing &&indexing, Filler &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert((bond.type() == "HOPUP") || (bond.type() == "HOPDN"));
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  idx_t size_up = indexing.size_ups();
  idx_t size_dn = indexing.size_dns();

  coeff_t traw = bond.coupling<coeff_t>();
  int s1 = bond[0];
  int s2 = bond[1];

  // Prepare bitmasks
  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  int l = std::min(s1, s2);
  int u = std::max(s1, s2);
  bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

  // Apply hoppings on upspins
  if (bond.type() == "HOPUP") {
    idx_t idx_up = 0;
    for (auto up : indexing.states_ups()) {
      if (bitops::popcnt(up & flipmask) == 1) {
        coeff_t t;
        if constexpr (is_complex<coeff_t>()) {
          t = (bitops::gbit(up, s1)) ? traw : hydra::conj(traw);
        } else {
          t = traw;
        }

        double fermi = bitops::popcnt(up & fermimask) & 1 ? -1. : 1.;
        coeff_t val = -t * fermi;

        bit_t up_flip = up ^ flipmask;
        idx_t idx_up_flip = indexing.index_ups(up_flip);
        for (idx_t idx_dn = 0; idx_dn < size_dn; ++idx_dn) {
          idx_t idx_out = idx_up_flip * size_dn + idx_dn;
          idx_t idx_in = idx_up * size_dn + idx_dn;
          fill(idx_out, idx_in, val);
        }
      }
      ++idx_up;
    }
  }

  // Apply hoppings on dnspins
  if (bond.type() == "HOPDN") {
    idx_t idx_dn = 0;
    for (auto dn : indexing.states_dns()) {
      if (bitops::popcnt(dn & flipmask) == 1) {
        coeff_t t;
        if constexpr (is_complex<coeff_t>()) {
          t = (bitops::gbit(dn, s1)) ? traw : hydra::conj(traw);
        } else {
          t = traw;
        }

        double fermi = bitops::popcnt(dn & fermimask) & 1 ? -1. : 1.;
        coeff_t val = -t * fermi;

        bit_t dn_flip = dn ^ flipmask;
        idx_t idx_dn_flip = indexing.index_dns(dn_flip);
        for (idx_t idx_up = 0; idx_up < size_up; ++idx_up) {
          idx_t idx_out = idx_up * size_dn + idx_dn_flip;
          idx_t idx_in = idx_up * size_dn + idx_dn;
          fill(idx_out, idx_in, val);
        }
      }
      ++idx_dn;
    }
  }
}

} // namespace hydra::electron
