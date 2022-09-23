#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void apply_hopping(Bond const &bond, Indexing &&indexing, Filler &&fill) {
  using combinatorics::Combinations;
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert((bond.type() == "HOPUP") || (bond.type() == "HOPDN"));
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  int n_sites = indexing.n_sites();
  int n_up = indexing.n_up();
  int n_holes = indexing.n_holes();
  int n_spins = indexing.n_spins();
  idx_t size_spins = indexing.size_spins();

  coeff_t traw = bond.coupling<coeff_t>();
  int s1 = bond[0];
  int s2 = bond[1];

  // Prepare bitmasks
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;
  bit_t spinsmask = ((bit_t)1 << n_spins) - 1;
  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  int l = std::min(s1, s2);
  int u = std::max(s1, s2);
  bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

  idx_t holes_idx = 0;
  for (auto holes : Combinations<bit_t>(n_sites, n_holes)) {

    // exactly one hole at sites s1, s2
    if (bitops::popcnt(holes & flipmask) & 1) {

      // get hopping parameter (complex conjugate if necessary)
      coeff_t t;
      if constexpr (is_complex<coeff_t>()) {
        t = (bitops::gbit(holes, s2)) ? traw : conj(traw);
      } else {
        t = traw;
      }

      bit_t not_holes = (~holes) & sitesmask;
      idx_t holes_offset = holes_idx * size_spins;
      bit_t new_holes = holes ^ flipmask;
      bit_t not_new_holes = (~new_holes) & sitesmask;
      idx_t new_holes_idx = indexing.index_holes(new_holes);
      idx_t new_holes_offset = new_holes_idx * size_spins;

      // Apply hoppings on upspins
      if (bond.type() == "HOPUP") {
        idx_t idx = holes_offset;
        for (auto spins : Combinations<bit_t>(n_spins, n_up)) {
          bit_t ups = bitops::deposit(spins, not_holes);
          if (bitops::popcnt(ups & flipmask) & 1) { // exactly one upspin
            bit_t new_ups = ups ^ flipmask;
            bit_t new_spins = bitops::extract(new_ups, not_new_holes);
            idx_t new_idx = new_holes_offset + indexing.index_spins(new_spins);
            bool fermi_sign = bitops::popcnt(ups & fermimask) & 1;

            // term -t \sum... (negative prefactor!)
            if (fermi_sign) {
              fill(new_idx, idx, t);
            } else {
              fill(new_idx, idx, -t);
            }
          }
          ++idx;
        }
      }

      // Apply hoppings on dnspins
      if (bond.type() == "HOPDN") {
        idx_t idx = holes_offset;
        for (auto spins : Combinations<bit_t>(n_spins, n_up)) {
          bit_t spins_neg = (~spins) & spinsmask;
          bit_t dns = bitops::deposit(spins_neg, not_holes);
          if (bitops::popcnt(dns & flipmask) & 1) { // exactly one dnspin
            bit_t new_dns = dns ^ flipmask;
            bit_t new_ups = ((~new_dns) & not_new_holes);
            bit_t new_spins = bitops::extract(new_ups, not_new_holes);
            idx_t new_idx = new_holes_offset + indexing.index_spins(new_spins);
            bool fermi_sign = bitops::popcnt(dns & fermimask) & 1;

            // term -t \sum... (negative prefactor!)
            if (fermi_sign) {
              fill(new_idx, idx, t);
            } else {
              fill(new_idx, idx, -t);
            }
          }
          ++idx;
        }
      } // if dn-hoppings
    }
    ++holes_idx;
  }
}

} // namespace hydra::tj
