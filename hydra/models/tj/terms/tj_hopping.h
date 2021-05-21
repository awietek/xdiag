#pragma once

#include <hydra/common.h>
#include <hydra/models/tj/tj.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>
#include <hydra/combinatorics/combinations.h>

namespace hydra::tjdetail {

template <class bit_t, class coeff_t, class Filler>
void do_hopping(BondList const &bonds, Couplings const &couplings,
                tJ<bit_t> const &block, Filler &&fill) {
  using utils::gbit;
  using utils::gbits;
  using utils::popcnt;

  int n_sites = block.n_sites();
  int n_up = block.n_up();

  auto hoppings = bonds.bonds_of_type("HOP") + bonds.bonds_of_type("TJHOP");
  auto hoppings_up =
      bonds.bonds_of_type("HOPUP") + bonds.bonds_of_type("TJHOPUP");
  auto hoppings_dn =
      bonds.bonds_of_type("HOPDN") + bonds.bonds_of_type("TJHOPDN");
  for (auto hop : hoppings + hoppings_up + hoppings_dn) {

    if (hop.size() != 2)
      HydraLog.err("Error computing tJ hopping: "
                   "hoppings must have exactly two sites defined");

    std::string cpl = hop.coupling();
    if (couplings.defined(cpl) && !lila::close(couplings[cpl], (complex)0.)) {

      int s1 = hop.site(0);
      int s2 = hop.site(1);
      bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
      int l = std::min(s1, s2);
      int u = std::max(s1, s2);
      bit_t spacemask = (((bit_t)1 << (u - l - 1)) - 1)
                        << (l + 1); // mask used for fermi sign

      // Apply hoppings on upspins
      if ((hop.type() == "HOP") || (hop.type() == "HOPUP") ||
          (hop.type() == "TJHOP") || (hop.type() == "TJHOPUP")) {
        idx_t idx_up = 0;
        for (auto up : Combinations<bit_t>(n_sites, n_up)) {
          if (popcnt(up & flipmask) == 1) {
            coeff_t t;
            if constexpr (is_complex<coeff_t>())
              t = (gbit(up, s1)) ? couplings[cpl] : lila::conj(couplings[cpl]);
            else
              t = lila::real(couplings[cpl]);

            double fermi = popcnt(up & spacemask) & 1 ? -1. : 1.;
            coeff_t val = -t * fermi;

            bit_t up_flip = up ^ flipmask;

            auto [dn_lower, dn_upper] = block.dn_limits_for_up_[idx_up];
            for (idx_t idx_in = dn_lower; idx_in < dn_upper; ++idx_in) {
              bit_t dn = block.dns_[idx_in];
              if ((dn & flipmask) == 0) { // no double occ possible
                idx_t idx_out = block.index(up_flip, dn);
                fill(idx_out, idx_in, val);
              }
            }
          }
          ++idx_up;
        }
      }

      // Apply hoppings on dnspins
      if ((hop.type() == "HOP") || (hop.type() == "HOPDN") ||
          (hop.type() == "TJHOP") || (hop.type() == "TJHOPDN")) {
        coeff_t t;
        if constexpr (is_complex<coeff_t>())
          t = couplings[cpl];
        else
          t = lila::real(couplings[cpl]);

        idx_t idx_up = 0;
        for (auto up : Combinations<bit_t>(n_sites, n_up)) {

          if ((up & flipmask) == 0) { // no double occ possible
            auto [dn_lower, dn_upper] = block.dn_limits_for_up_[idx_up];
            for (idx_t idx_in = dn_lower; idx_in < dn_upper; ++idx_in) {
              bit_t dn = block.dns_[idx_in];

              if (popcnt(dn & flipmask) == 1) {
                bit_t dn_flip = dn ^ flipmask;
                auto it =
                    std::lower_bound(block.dns_.begin() + dn_lower,
                                     block.dns_.begin() + dn_upper, dn_flip);
                idx_t idx_out = std::distance(block.dns_.begin(), it);

                // Complex conjugate for complex coefficients
                if constexpr (is_complex<coeff_t>()) {
                  if (gbit(dn, s2)) {
                    // Take fermi sign into account
                    if (popcnt(dn & spacemask) & 1) {
                      fill(idx_out, idx_in, lila::conj(t));
                    } else {
                      fill(idx_out, idx_in, -lila::conj(t));
                    }
                  } else {
                    if (popcnt(dn & spacemask) & 1) {
                      fill(idx_out, idx_in, t);
                    } else {
                      fill(idx_out, idx_in, -t);
                    }
                  }
                  // Real case doesnt need complex conjugate
                } else {
                  if (popcnt(dn & spacemask) & 1) {
                    fill(idx_out, idx_in, t);
                  } else {
                    fill(idx_out, idx_in, -t);
                  }
                }
              }
            }
          }
	  ++idx_up;
        }
      }
    }
  }
}
} // namespace hydra::tjdetail
