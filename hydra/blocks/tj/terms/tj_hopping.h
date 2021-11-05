#pragma once

#include <lila/utils/logger.h>

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/tj/tj_utils.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms::tj {

template <class bit_t, class coeff_t, class Filler>
void do_hopping(BondList const &bonds, Couplings const &couplings,
                tJ<bit_t> const &block, Filler &&fill) {
  using bitops::bits_to_string;
  using bitops::deposit;
  using bitops::extract;
  using bitops::gbit;
  using bitops::popcnt;
  using combinatorics::Combinations;

  int n_sites = block.n_sites();
  int N = n_sites;
  int nup = block.n_up();
  int ndn = block.n_dn();
  int n_holes = n_sites - nup - ndn;
  int charge = nup + ndn;
  int size_spins = block.size_spins_;

  auto hoppings = bonds.bonds_of_type("HOP") + bonds.bonds_of_type("TJHOP");
  auto hoppings_up =
      bonds.bonds_of_type("HOPUP") + bonds.bonds_of_type("TJHOPUP");
  auto hoppings_dn =
      bonds.bonds_of_type("HOPDN") + bonds.bonds_of_type("TJHOPDN");
  for (auto hop : hoppings + hoppings_up + hoppings_dn) {

    if (hop.size() != 2)
      lila::Log.err("Error computing tJ hopping: "
                    "hoppings must have exactly two sites defined");

    std::string cpl = hop.coupling();
    if (couplings.defined(cpl) && !lila::close(couplings[cpl], (complex)0.)) {

      bit_t sitesmask = ((bit_t)1 << n_sites) - 1;
      int s1 = hop[0];
      int s2 = hop[1];
      bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
      int l = std::min(s1, s2);
      int u = std::max(s1, s2);
      bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

      idx_t holes_idx = 0;
      for (auto holes : Combinations<bit_t>(n_sites, n_holes)) {

        if (popcnt(holes & flipmask) & 1) { // exactly one hole at sites s1, s2

          // get hopping parameter (complex conjugate if necessary)
          coeff_t t;
          if constexpr (is_complex<coeff_t>()) {
            t = (gbit(holes, s2)) ? couplings[cpl] : lila::conj(couplings[cpl]);
          } else {
            t = lila::real(couplings[cpl]);
          }

          bit_t not_holes = (~holes) & sitesmask;
          idx_t holes_offset = holes_idx * size_spins;
          bit_t new_holes = holes ^ flipmask;
	  bit_t not_new_holes = (~new_holes) & sitesmask;
          idx_t new_holes_idx = block.lintable_holes_.index(new_holes);
          idx_t new_holes_offset = new_holes_idx * size_spins;

          // Apply hoppings on upspins
          if ((hop.type() == "HOP") || (hop.type() == "HOPUP") ||
              (hop.type() == "TJHOP") || (hop.type() == "TJHOPUP")) {
            idx_t idx = holes_offset;
            for (auto spins : Combinations<bit_t>(charge, nup)) {
              bit_t ups = bitops::deposit(spins, not_holes);
              if (popcnt(ups & flipmask) & 1) { // exactly one upspin
                bit_t new_ups = ups ^ flipmask;
                bit_t new_spins = bitops::extract(new_ups, not_new_holes);
                idx_t new_idx =
                    new_holes_offset + block.lintable_spins_.index(new_spins);
                bool fermi_sign = popcnt(ups & fermimask) & 1;

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
          if ((hop.type() == "HOP") || (hop.type() == "HOPDN") ||
              (hop.type() == "TJHOP") || (hop.type() == "TJHOPDN")) {
            idx_t idx = holes_offset;
            for (auto spins : Combinations<bit_t>(charge, nup)) {
              bit_t spins_neg = (~spins) & sitesmask;
              bit_t dns = bitops::deposit(spins_neg, not_holes);
              if (popcnt(dns & flipmask) & 1) { // exactly one dnspin
                bit_t new_dns = dns ^ flipmask;
                bit_t new_ups = (~new_dns) & sitesmask;
                bit_t new_spins = bitops::extract(new_ups, not_new_holes);
                idx_t new_idx =
                    new_holes_offset + block.lintable_spins_.index(new_spins);
                bool fermi_sign = popcnt(dns & fermimask) & 1;

		// lila::Log("holes : {} spins : {}", bits_to_string(holes, N), bits_to_string(spins, charge));
		// lila::Log("mask  : {}", bits_to_string(flipmask, N));
		// lila::Log("nholes: {}", bits_to_string(new_holes, N));
		// lila::Log("dns   : {}", bits_to_string(dns, N));
		// lila::Log("mask  : {}", bits_to_string(flipmask, N));
		// lila::Log("ndns  : {}", bits_to_string(new_dns, N));
		// lila::Log("nholes: {} nspins: {}\n", bits_to_string(new_holes, N), bits_to_string(new_spins, charge));

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
  }
}
} // namespace hydra::terms::tj
