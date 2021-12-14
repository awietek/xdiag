#pragma once

#include <hydra/common.h>

#include <hydra/blocks/utils/block_utils.h>

#include <hydra/bitops/bitops.h>

#include <hydra/indexing/electron/electron_symmetric_indexing.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void electron_symmetric_hopping(BondList const &bonds,
                                Couplings const &couplings, Indexing &&indexing,
                                Filler &&fill) {
  using bitops::gbit;
  using bitops::popcnt;

  // Get group/irrep info
  auto const &group_action = indexing.group_action();
  auto const &irrep = indexing.irrep();
  std::vector<coeff_t> bloch_factors;
  if constexpr (is_complex<coeff_t>()) {
    bloch_factors = irrep.characters();
  } else {
    bloch_factors = irrep.characters_real();
  }
  std::vector<coeff_t> prefacs(irrep.size());

  assert(group_action.n_symmetries() == bloch_factors.size());

  auto clean_bonds =
      utils::clean_bondlist(bonds, couplings, {"HOP", "HOPUP", "HOPDN"}, 2);
  for (auto bond : clean_bonds) {

    std::string type = bond.type();
    std::string cpl = bond.coupling();

    utils::check_sites_disjoint(bond);
    int s1 = bond[0];
    int s2 = bond[1];

    auto [t, tconj] = utils::get_coupling_and_conj<coeff_t>(couplings, cpl);

    // Prepare bitmasks
    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    int l = std::min(s1, s2);
    int u = std::max(s1, s2);
    bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

    coeff_t val = 0;

    // Apply hoppings on dnspins
    if ((type == "HOP") || (type == "HOPDN")) {

      for (idx_t idx_up = 0; idx_up < indexing.n_rep_ups(); ++idx_up) {
        bit_t ups = indexing.rep_ups(idx_up);
        idx_t up_offset = indexing.ups_offset(idx_up);
        auto syms = indexing.syms_ups(ups);
        auto dnss = indexing.dns_for_ups_rep(ups);
        auto norms = indexing.norms_for_ups_rep(ups);

        // trivial up-stabilizer (likely)
        if (syms.size() == 1) {
          idx_t idx_dn = 0;
          for (bit_t dns : dnss) {

            if (popcnt(dns & flipmask) == 1) {

              idx_t idx_in = up_offset + idx_dn;
              bit_t dns_flip = dns ^ flipmask;
              idx_t idx_dn_flip = indexing.index_dns(dns_flip);
              idx_t idx_out = up_offset + idx_dn_flip;
              bool fermi = popcnt(dns & fermimask) & 1;

              // Complex conjugate t if necessary
              if constexpr (is_complex<coeff_t>()) {
                val = -(gbit(dns, s1) ? t : tconj);
              } else {
                val = -t;
              } // Comment: norm is always 1.0 for trivial stabilizers

              fill(idx_out, idx_in, fermi ? -val : val);
            }

            ++idx_dn;
          }

          // non-trivial up-stabilizer (unlikely)
        } else {
          idx_t idx_dns = 0;
          for (bit_t dns : dnss) {

            if (popcnt(dns & flipmask) == 1) {

              idx_t idx_in = up_offset + idx_dns;
              bit_t dns_flip = dns ^ flipmask;
              auto [idx_dns_flip, fermi_dn, sym] =
                  indexing.index_dns_fermi_sym(dns_flip, syms, dnss, fermimask);

              if (idx_dns_flip != invalid_index) {
                idx_t idx_out = up_offset + idx_dns_flip;
                bool fermi_up = indexing.fermi_bool_ups(sym, ups);

                // Complex conjugate t if necessary
                if constexpr (is_complex<coeff_t>()) {
                  val = -(gbit(dns, s1) ? t : tconj) * bloch_factors[sym] *
                        norms[idx_dns_flip] / norms[idx_dns];
                } else {
                  val = -t * bloch_factors[sym] * norms[idx_dns_flip] /
                        norms[idx_dns];
                }

                fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
              }
            }

            ++idx_dns;
          }
        } // if trivial stabilizer or not
      }   // loop over ups
    }     // if dn-hoppings

    // Apply hoppings on upspins
    if ((type == "HOP") || (type == "HOPUP")) {
      for (idx_t idx_up = 0; idx_up < indexing.n_rep_ups(); ++idx_up) {
        bit_t ups = indexing.rep_ups(idx_up);

        if (popcnt(ups & flipmask) != 1)
          continue;

        bit_t ups_flip = ups ^ flipmask;
        idx_t idx_ups_flip = indexing.index_ups(ups_flip);
        bit_t ups_flip_rep = indexing.rep_ups(idx_ups_flip);

        // Get limits, syms, and dns for ingoing ups
        idx_t ups_offset_in = indexing.ups_offset(idx_up);
        auto dnss_in = indexing.dns_for_ups_rep(ups);
        auto norms_in = indexing.norms_for_ups_rep(ups);

        // Get limits, syms, and dns for outgoing ups
        idx_t ups_offset_out = indexing.ups_offset(idx_ups_flip);
        auto syms_ups_out = indexing.syms_ups(ups_flip);
        auto dnss_out = indexing.dns_for_ups_rep(ups_flip_rep);
        auto norms_out = indexing.norms_for_ups_rep(ups_flip_rep);

        // trivial up-stabilizer (likely)
        if (syms_ups_out.size() == 1) {
          int sym = syms_ups_out.front();

          // Complex conjugate t if necessary
          coeff_t prefac;
          if constexpr (is_complex<coeff_t>()) {
            prefac = -(gbit(ups, s1) ? t : tconj) * bloch_factors[sym];
          } else {
            prefac = -t * bloch_factors[sym];
          }

          // Fermi-sign of up spins
          bool fermi_up = (popcnt(ups & fermimask) & 1);
          fermi_up ^= indexing.fermi_bool_ups(sym, ups_flip);

          idx_t idx_dn = 0;
          for (bit_t dns : dnss_in) {
            idx_t idx_in = ups_offset_in + idx_dn;
            bit_t dns_rep = group_action.apply(sym, dns);
            idx_t idx_out = ups_offset_out + indexing.index_dns(dns_rep);
            bool fermi_dn = indexing.fermi_bool_dns(sym, dns);

            val = prefac / norms_in[idx_dn]; // norms_out = 1.0 in this case

            fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);

            ++idx_dn;
          }
          // non-trivial up-stabilizer (unlikely)
        } else {
          bool fermi_up_hop = (popcnt(ups & fermimask) & 1);
          auto syms = syms_ups_out;

          // Fix the bloch/prefactors
          if constexpr (is_complex<coeff_t>()) {
            for (int i = 0; i < (int)irrep.size(); ++i) {
              prefacs[i] = -(gbit(ups, s1) ? t : tconj) * irrep.character(i);
            }
          } else {
            for (int i = 0; i < (int)irrep.size(); ++i) {
              prefacs[i] = -t * lila::real(irrep.character(i));
            }
          }

          idx_t idx_dn = 0;
          for (bit_t dns : dnss_in) {
            idx_t idx_in = ups_offset_in + idx_dn;

            auto [idx_dn_out, fermi_dn, sym] =
                indexing.index_dns_fermi_sym(dns, syms, dnss_out);

            if (idx_dn_out != invalid_index) {
              idx_t idx_out = ups_offset_out + idx_dn_out;
              bool fermi_up =
                  fermi_up_hop ^ indexing.fermi_bool_ups(sym, ups_flip);
              coeff_t val =
                  prefacs[sym] * norms_out[idx_dn_out] / norms_in[idx_dn];

              fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
            }
            ++idx_dn;
          }
        } // if trivial stabilizer or not
      }   // loop over ups
    }     // if up-hoppings
  }       // loop over bonds
}
} // namespace hydra::terms
