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
  using bitops::gbit;
  using bitops::popcnt;
  using combinatorics::Combinations;

  int n_sites = block.n_sites();
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

      int s1 = hop[0];
      int s2 = hop[1];
      int l = std::min(s1, s2);
      bit_t holemask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

      idx_t holes_idx = 0;
      for (auto holes : Combinations<bit_t>(n_sites, n_holes)) {

        if (popcnt(holes & holemask) == 1) {

          // get hopping parameter (complex conjugate if necessary)
          coeff_t t;
          if constexpr (is_complex<coeff_t>()) {
            t = (gbit(holes, s2)) ? couplings[cpl] : lila::conj(couplings[cpl]);
          } else {
            t = lila::real(couplings[cpl]);
          }

          bit_t new_holes = holes ^ holemask;
          bit_t new_holes_idx = block.lintable_holes_.index(new_holes);

          idx_t holes_offset = holes_idx * size_spins;
          idx_t new_holes_offset = new_holes_idx * size_spins;

          // new site number for spins
          int s1c = utils::tj_site_compressed(holes, s1);
          int s2c = utils::tj_site_compressed(holes, s2);
          int lc = std::min(s1c, s2c);
          int uc = std::max(s1c, s2c);

          // constant masks for bit manipulation
          bit_t sitesmask = (((bit_t)1 << (uc - lc + 1)) - 1) << lc;
          bit_t sitesmask_neg = ~sitesmask;
          bit_t fermimask = (((bit_t)1 << (uc - lc - 1)) - 1) << (lc + 1);
          bit_t up_at_lc = (bit_t)1 << lc;
          bit_t up_at_uc = (bit_t)1 << uc;

          // Apply hoppings on upspins
          if ((hop.type() == "HOP") || (hop.type() == "HOPUP") ||
              (hop.type() == "TJHOP") || (hop.type() == "TJHOPUP")) {
            if (gbit(holes, l)) { // up-electron hopping u -> l
              idx_t idx = holes_offset;
              for (auto spins : Combinations<bit_t>(charge, nup)) {
                if (spins & up_at_uc) {
                  bit_t spins_canvas = spins & sitesmask_neg;
                  bit_t spins_center = spins & sitesmask;
                  bit_t spins_center_rot =
                      (spins_center << 1 | up_at_lc) & sitesmask;
                  bit_t new_spins = spins_canvas & spins_center_rot;
                  bit_t new_spins_idx = block.lintable_spins_.index(new_spins);
                  bit_t new_idx = new_holes_offset + new_spins_idx;
                  if (popcnt(spins & fermimask) & 1) { // fermi sign
                    fill(new_idx, idx, t);
                  } else {
                    fill(new_idx, idx, -t);
                  }
                }
                ++idx;
              }
            } else { // up-electron hopping l -> u
              idx_t idx = holes_offset;
              for (auto spins : Combinations<bit_t>(charge, nup)) {
                if (spins & up_at_lc) {
                  bit_t spins_canvas = spins & sitesmask_neg;
                  bit_t spins_center = spins & sitesmask;
                  bit_t spins_center_rot =
                      (spins_center >> 1 | up_at_uc) & sitesmask;
                  bit_t new_spins = spins_canvas & spins_center_rot;
                  bit_t new_spins_idx = block.lintable_spins_.index(new_spins);
                  bit_t new_idx = new_holes_offset + new_spins_idx;
                  if (popcnt(spins & fermimask) & 1) { // fermi sign
                    fill(new_idx, idx, t);
                  } else {
                    fill(new_idx, idx, -t);
                  }
                }
                ++idx;
              }
            }
          } // if up-hoppings

          if ((hop.type() == "HOP") || (hop.type() == "HOPDN") ||
              (hop.type() == "TJHOP") || (hop.type() == "TJHOPDN")) {
            if (gbit(holes, l)) { // dn-electron hopping u -> l
              idx_t idx = holes_offset;
              for (auto spins : Combinations<bit_t>(charge, nup)) {
                if ((~spins) & up_at_uc) {
                  bit_t spins_canvas = spins & sitesmask_neg;
                  bit_t spins_center = spins & sitesmask;
                  bit_t spins_center_rot = (spins_center << 1) & sitesmask;
                  bit_t new_spins = spins_canvas & spins_center_rot;
                  bit_t new_spins_idx = block.lintable_spins_.index(new_spins);
                  bit_t new_idx = new_holes_offset + new_spins_idx;
                  if (popcnt((~spins) & fermimask) & 1) { // fermi sign
                    fill(new_idx, idx, t);
                  } else {
                    fill(new_idx, idx, -t);
                  }
                }
                ++idx;
              }
            } else { // up-electron hopping l -> u
              idx_t idx = holes_offset;
              for (auto spins : Combinations<bit_t>(charge, nup)) {
                if ((~spins) & up_at_lc) {
                  bit_t spins_canvas = spins & sitesmask_neg;
                  bit_t spins_center = spins & sitesmask;
                  bit_t spins_center_rot = (spins_center >> 1) & sitesmask;
                  bit_t new_spins = spins_canvas & spins_center_rot;
                  bit_t new_spins_idx = block.lintable_spins_.index(new_spins);
                  bit_t new_idx = new_holes_offset + new_spins_idx;
                  if (popcnt((~spins) & fermimask) & 1) { // fermi sign
                    fill(new_idx, idx, t);
                  } else {
                    fill(new_idx, idx, -t);
                  }
                }
                ++idx;
              }
            }
          } // if dn-hoppings
        }
        ++holes_idx;
      }
    }
  }
}
} // namespace hydra::terms::tj

// int s1c = utils::tj_site_compressed(holes, s1);
// int s2c = utils::tj_site_compressed(holes, s2);

//   // Apply hoppings on upspins
//   if ((hop.type() == "HOP") || (hop.type() == "HOPUP") ||
//       (hop.type() == "TJHOP") || (hop.type() == "TJHOPUP")) {

//     idx_t idx_up = 0;
//     for (auto up : Combinations<bit_t>(n_sites, nup)) {
//       if (popcnt(up & flipmask) == 1) {

//         // Complex conjugate if necessary
//         coeff_t t;
//         if constexpr (is_complex<coeff_t>())
//           t = (gbit(up, s1)) ? couplings[cpl] :
//           lila::conj(couplings[cpl]);
//         else
//           t = lila::real(couplings[cpl]);

//         // Compute the fermi sign
//         double fermi = popcnt(up & spacemask) & 1 ? -1. : 1.;
//         coeff_t val = -t * fermi;

//         // Compute flipped ups and new up idx
//         bit_t up_flip = up ^ flipmask;
//         idx_t idx_up_flip = block.lintable_up_.index(up_flip) * size_dn;

//         // Define compressed mask for dn spins
//         int s1c = utils::tj_site_compressed(up, s1);
//         int s2c = utils::tj_site_compressed(up, s2);
//         bit_t cflipmask = ((bit_t)1 << s1c) | ((bit_t)1 << s2c);

//         // Go over all dn configurations and fill
//         for (idx_t idx_dn = 0; idx_dn < size_dn; ++idx_in) {
//           bit_t dnc = block.dncs_[idx_in];
//           if ((dn & cflipmask) == 0) { // no double occ possible
//             idx_t idx_in = idx_up + idx_dn;
//             idx_t idx_out = idx_up_flip + idx_dn;
//             fill(idx_out, idx_in, val);
//           }
//         }
//       }
//       idx_up += size_dn;
//     }
//   }

//   // Apply hoppings on dnspins
//   if ((hop.type() == "HOP") || (hop.type() == "HOPDN") ||
//       (hop.type() == "TJHOP") || (hop.type() == "TJHOPDN")) {
//     coeff_t t;
//     if constexpr (is_complex<coeff_t>())
//       t = couplings[cpl];
//     else
//       t = lila::real(couplings[cpl]);

//     idx_t idx_dn = 0;
//     for (auto dn : Combinations<bit_t>(n_sites, ndn)) {
//       if (popcnt(dn & flipmask) == 1) {

//         // Complex conjugate if necessary
//         coeff_t t;
//         if constexpr (is_complex<coeff_t>())
//           t = (gbit(dn, s1)) ? couplings[cpl] :
//           lila::conj(couplings[cpl]);
//         else
//           t = lila::real(couplings[cpl]);

//         // Compute the fermi sign
//         double fermi = popcnt(dn & spacemask) & 1 ? -1. : 1.;
//         coeff_t val = -t * fermi;

//         // Compute flipped dns and new up idx
//         bit_t dn_flip = dn ^ flipmask;

//         // Define compressed mask for up spins
//         int s1c = utils::tj_site_compressed(up, s1);
//         int s2c = utils::tj_site_compressed(up, s2);
//         bit_t cflipmask = ((bit_t)1 << s1c) | ((bit_t)1 << s2c);

// 	    // Go over all up configurations and fill
//         for (idx_t idx_up = 0; idx_up < size_up; ++idx_up) {
//           bit_t upc = block.upcs_[idx_up];
//           if ((upc & cflipmask) == 0) { // no double occ possible
// 		auto [up, dnc] = utils::upc_dn_to_up_dnc(upc, dn);
// 		idx_t up_offset = lintable_up_.index(up) * size_dn;
//             idx_t idx_in = up_offset + idx_dn;
// 		idx_t idx_dn_flip = lintable_dn_.index
//             idx_t idx_out = idx_dn;
//             fill(idx_out, idx_in, val);
//           }
// 	      up_offset += size_dn;
//         }              bit_t dn = block.dns_[idx_in];

//       }
// 	  ++idx_dn;
//     }
//   }
// }
// if ((dn & flipmask) == 0) { // no double occ possible
//   auto [dn_lower, dn_upper] = block.dn_limits_for_up_[idx_up];
//   for (idx_t idx_in = dn_lower; idx_in < dn_upper; ++idx_in) {
//     bit_t dn = block.dns_[idx_in];

//     if (popcnt(dn & flipmask) == 1) {
//       bit_t dn_flip = dn ^ flipmask;
//       auto it =
//           std::lower_bound(block.dns_.begin() + dn_lower,
//                            block.dns_.begin() + dn_upper, dn_flip);
//       idx_t idx_out = std::distance(block.dns_.begin(), it);

//       // Complex conjugate for complex coefficients
//       if constexpr (is_complex<coeff_t>()) {
//         if (gbit(dn, s2)) {
//           // Take fermi sign into account
//           if (popcnt(dn & spacemask) & 1) {
//             fill(idx_out, idx_in, lila::conj(t));
//           } else {
//             fill(idx_out, idx_in, -lila::conj(t));
//           }
//         } else {
//           if (popcnt(dn & spacemask) & 1) {
//             fill(idx_out, idx_in, t);
//           } else {
//             fill(idx_out, idx_in, -t);
//           }
//         }
//         // Real case doesnt need complex conjugate
//       } else {
//         if (popcnt(dn & spacemask) & 1) {
//           fill(idx_out, idx_in, t);
//         } else {
//           fill(idx_out, idx_in, -t);
//         }
//       }
//     }
//   }
// }
// ++idx_up;
