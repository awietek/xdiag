#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

#include <hydra/blocks/electron/terms/compile.h>

#include <hydra/blocks/electron/terms/apply_exchange_sym.h>
#include <hydra/blocks/electron/terms/apply_hopping_sym.h>
#include <hydra/blocks/electron/terms/apply_ising_sym.h>
#include <hydra/blocks/electron/terms/apply_number_sym.h>
#include <hydra/blocks/electron/terms/apply_u_sym.h>

#include <hydra/blocks/electron/terms/apply_exchange_no_sym.h>
#include <hydra/blocks/electron/terms/apply_hopping_no_sym.h>
#include <hydra/blocks/electron/terms/apply_ising_no_sym.h>
#include <hydra/blocks/electron/terms/apply_number_no_sym.h>
#include <hydra/blocks/electron/terms/apply_u_no_sym.h>

#include <hydra/blocks/electron/terms/apply_generic_term_diag_no_sym.h>
#include <hydra/blocks/electron/terms/apply_generic_term_dns_no_sym.h>
#include <hydra/blocks/electron/terms/apply_generic_term_ups_no_sym.h>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class Fill>
void apply_terms(BondList const &bonds, IndexingIn const &indexing_in,
                 IndexingOut const &indexing_out, Fill &fill) {

  auto bonds_compiled = electron::compile(bonds);

  for (Bond bond : bonds_compiled) {
    std::string type = bond.type();

    // Hopping terms
    if ((type == "HOPUP") || (type == "HOPDN")) {
      if (bond.size() != 2) {
        Log.err("Error trying to apply hopping: number of sites found to "
                "be {} instead of 2",
                bond.size());
      }
      int s1 = bond[0];
      int s2 = bond[1];
      if (s1 == s2) {
        Log.err("Error trying to apply hopping: site 1 and site 2 must be "
                "distinct from one another. If intentional, use the number "
                "operator instead");
      }
      bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
      int l = std::min(s1, s2);
      int u = std::max(s1, s2);
      bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);
      coeff_t t = bond.coupling<coeff_t>();

      auto non_zero_term = [&flipmask](bit_t const &spins) -> bool {
        return bitops::popcnt(spins & flipmask) & 1;
      };
      auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
        bool fermi = bitops::popcnt(spins & fermimask) & 1;
        spins ^= flipmask;
        if constexpr (is_complex<coeff_t>()) {
          coeff_t tt = (bitops::gbit(spins, s1)) ? t : hydra::conj(t);
          return {spins, fermi ? tt : -tt};
        } else {
          return {spins, fermi ? t : -t};
        }
      };

      if (type == "HOPUP") {
        if constexpr (symmetric) {
          electron::apply_hopping_sym<bit_t, coeff_t>(bond, indexing_in, fill);
        } else {
          electron::apply_generic_term_ups_no_sym<bit_t, coeff_t>(
              indexing_in, indexing_out, non_zero_term, term_action, fill);
        }
      } else if (type == "HOPDN") {
        if constexpr (symmetric) {
          electron::apply_hopping_sym<bit_t, coeff_t>(bond, indexing_in, fill);
        } else {
          electron::apply_generic_term_dns_no_sym<bit_t, coeff_t, false>(
              indexing_in, indexing_out, non_zero_term, term_action, fill);
        }
      }

    } // if ((type == "HOPUP") || (type == "HOPDN"))

    else if (type == "EXCHANGE") {
      if constexpr (symmetric) {
        electron::apply_exchange_sym<bit_t, coeff_t>(bond, indexing_in, fill);
      } else {
        electron::apply_exchange_no_sym<bit_t, coeff_t>(bond, indexing_in,
                                                        fill);
      }
    } // if (type == "EXCHANGE")

    else if (type == "ISING") {
      if constexpr (symmetric) {
        electron::apply_ising_sym<bit_t, coeff_t>(bond, indexing_in, fill);
      } else {
        electron::apply_ising_no_sym<bit_t, coeff_t>(bond, indexing_in, fill);
      }
    } // if (type == "ISING")

    else if ((type == "NUMBERUP") || (type == "NUMBERDN")) {
      if constexpr (symmetric) {
        electron::apply_number_sym<bit_t, coeff_t>(bond, indexing_in, fill);
      } else {
        electron::apply_number_no_sym<bit_t, coeff_t>(bond, indexing_in, fill);
      }
    } //  if ((type == "NUMBERUP") || (type == "NUMBERDN"))

    // Raising / lowering operators
    else if ((type == "CDAGUP") || (type == "CDAGDN") || (type == "CUP") ||
             (type == "CDN")) {
      if (bond.size() != 1) {
        Log.err("Error trying to apply Cdagup/Cdagdn/Cup/Cdn: number of sites "
                "found to be {} instead of 1",
                bond.size());
      }
      int s = bond[0];
      bit_t site_mask = (bit_t)1 << s;
      bit_t fermi_mask = site_mask - 1;
      coeff_t c = bond.coupling<coeff_t>();

      // Raising operators
      if ((type == "CDAGUP") || (type == "CDAGDN")) {
        auto non_zero_term = [&](bit_t const &spins) -> bool {
          return (spins & site_mask) == 0;
        };
        auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
          bool fermi = bitops::popcnt(spins & fermi_mask) & 1;
          return {spins ^ site_mask, fermi ? -c : c};
        };

        if (type == "CDAGUP") {
          if constexpr (symmetric) {
          } else {
            electron::apply_generic_term_ups_no_sym<bit_t, coeff_t>(
                indexing_in, indexing_out, non_zero_term, term_action, fill);
          }
        } else if (type == "CDAGDN") {
          if constexpr (symmetric) {
          } else {
            electron::apply_generic_term_dns_no_sym<bit_t, coeff_t, true>(
                indexing_in, indexing_out, non_zero_term, term_action, fill);
          }
        }

        // Lowering operators
      } else if ((type == "CUP") || (type == "CDN")) {
        auto non_zero_term = [&](bit_t const &spins) -> bool {
          return (spins & site_mask);
        };
        auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
          bool fermi = bitops::popcnt(spins & fermi_mask) & 1;
	  return {spins ^ site_mask, fermi ? -c : c};
        };

        if (type == "CUP") {
          if constexpr (symmetric) {
          } else {
            electron::apply_generic_term_ups_no_sym<bit_t, coeff_t>(
                indexing_in, indexing_out, non_zero_term, term_action, fill);
          }
        } else if (type == "CDN") {
          if constexpr (symmetric) {
          } else {
            electron::apply_generic_term_dns_no_sym<bit_t, coeff_t, true>(
                indexing_in, indexing_out, non_zero_term, term_action, fill);
          }
        }
      }
    }

    else {
      Log.err("Error in electron::apply_terms: Unknown bond type {}", type);
    }

  } // for (Bond bond : bonds)

  if (bonds.coupling_defined("U")) {
    coeff_t U = bonds.coupling<coeff_t>("U");
    // auto term_action = [&U](bit_t const &ups, bit_t const &dns) {
    //   return U * (double)bitops::popcnt(ups & dns);
    // };

    if constexpr (symmetric) {
      electron::apply_u_sym<bit_t, coeff_t>(U, indexing_in, fill);
    } else {
      electron::apply_u_no_sym<bit_t, coeff_t>(U, indexing_in, fill);
    }
  }
}

} // namespace hydra::electron
