#pragma once

#include <lila/utils/logger.h>

#include <hydra/common.h>
#include <hydra/models/electron_symmetric/electron_symmetric.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/symmetry_utils.h>

#include <hydra/models/utils/model_utils.h>
#include <hydra/utils/bitops.h>

namespace hydra::electronterms {

template <class bit_t, class coeff_t, class GroupAction, class Filler>
void do_hopping_symmetric(BondList const &bonds, Couplings const &couplings,
                          ElectronSymmetric<bit_t, GroupAction> const &block,
                          Filler &&fill) {
  using utils::gbit;
  using utils::popcnt;

  auto const &group_action = block.group_action();
  auto const &irrep = block.irrep();

  auto hoppings = bonds.bonds_of_type("HOP");
  auto hoppings_up = bonds.bonds_of_type("HOPUP");
  auto hoppings_dn = bonds.bonds_of_type("HOPDN");
  for (auto hop : hoppings + hoppings_up + hoppings_dn) {

    if (hop.size() != 2)
      lila::Log.err("Error computing Electron hopping: "
                    "hoppings must have exactly two sites defined");

    std::string cpl = hop.coupling();
    if (utils::coupling_is_non_zero(hop, couplings)) {
      int s1 = hop.site(0);
      int s2 = hop.site(1);
      bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
      int l = std::min(s1, s2);
      int u = std::max(s1, s2);
      bit_t spacemask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

      // Define hopping amplitude
      coeff_t t = 0;
      coeff_t tconj = 0;
      if constexpr (is_complex<coeff_t>()) {
        t = couplings[cpl];
        tconj = lila::conj(t);
      } else {
        t = lila::real(couplings[cpl]);
        tconj = lila::real(couplings[cpl]);
      }
      coeff_t val = 0;

      // Apply hoppings on dnspins
      if ((hop.type() == "HOP") || (hop.type() == "HOPDN")) {

        idx_t idx_up = 0;
        for (bit_t ups : block.reps_up_) {
          auto const &dnss = block.dns_for_up_rep(ups);
          auto const &norms = block.norms_for_up_rep(ups);

          idx_t up_offset = block.up_offsets_[idx_up];
          auto [sym_lower, sym_upper] = block.sym_limits_up(ups);
          assert(sym_upper - sym_lower > 0);

          // trivial up-stabilizer (likely)
          if (sym_upper - sym_lower == 1) {
            idx_t idx_dn = 0;
            for (bit_t dns : dnss) {

              // If hopping is possible
              if (popcnt(dns & flipmask) == 1) {
                idx_t idx_in = up_offset + idx_dn;
                bit_t dns_flip = dns ^ flipmask;
                idx_t idx_dn_flip = block.lintable_dns_.index(dns_flip);
                idx_t idx_out = up_offset + idx_dn_flip;
                double fermi_hop = popcnt(dns & spacemask) & 1 ? -1. : 1.;

                // Complex conjugate t if necessary
                if constexpr (is_complex<coeff_t>()) {
                  val = -(gbit(dns, s1) ? t : tconj) * fermi_hop *
                        norms[idx_dn_flip] / norms[idx_dn];
                } else {
                  val = -t * fermi_hop * norms[idx_dn_flip] / norms[idx_dn];
                }

                fill(idx_out, idx_in, val);
              }
              ++idx_dn;
            }
            // non-trivial up-stabilizer (unlikely)
          } else {
            idx_t idx_dn = 0;
            for (bit_t dns : dnss) {

              // If hopping is possible
              if (popcnt(dns & flipmask) == 1) {
                idx_t idx_in = up_offset + idx_dn;
                bit_t dns_flip = dns ^ flipmask;

                // Determine dns representative
                bit_t rep_dns = std::numeric_limits<bit_t>::max();
                int rep_sym = 0;
                for (idx_t sym_idx = sym_lower; sym_idx < sym_upper;
                     ++sym_idx) {
                  int sym = block.syms_up_[sym_idx];
                  bit_t tdns = group_action.apply(sym, dns_flip);
                  if (tdns < rep_dns) {
                    rep_dns = tdns;
                    rep_sym = sym;
                  }
                }

                // Find index of rep dns
                auto it = std::lower_bound(dnss.begin(), dnss.end(), rep_dns);
                if ((it != dnss.end()) && (*it == rep_dns)) {
                  idx_t idx_dn_flip = std::distance(dnss.begin(), it);
                  idx_t idx_out = up_offset + idx_dn_flip;
                  double fermi_hop = popcnt(dns & spacemask) & 1 ? -1. : 1.;
                  if (block.fermi_bool_up(rep_sym, ups) !=
                      block.fermi_bool_dn(rep_sym, dns_flip))
                    fermi_hop = -fermi_hop;

                  // Complex conjugate t if necessary
                  if constexpr (is_complex<coeff_t>()) {
                    val = -(gbit(dns, s1) ? t : tconj) * fermi_hop *
                          irrep.character(rep_sym) * norms[idx_dn_flip] /
                          norms[idx_dn];
                  } else {
                    val = -t * fermi_hop *
                          lila::real(irrep.character(rep_sym)) *
                          norms[idx_dn_flip] / norms[idx_dn];
                  }

                  // // Debug
                  // idx_t idx_out2 = block.index(ups, rep_dns);
                  // assert(idx_out == idx_out2);
                  // int n_sites = block.n_sites();
                  // lila::Log.out("s1: {} s2: {} {};{} -> {};{} -> {};{}
                  // idx_in: {}, idx_out: {}, val: {:.4f}", s1, s2,
                  //               bits_to_string(ups, n_sites),
                  //               bits_to_string(dns, n_sites),
                  //               bits_to_string(ups, n_sites),
                  //               bits_to_string(dns_flip, n_sites),
                  //               bits_to_string(ups, n_sites),
                  //               bits_to_string(rep_dns, n_sites), idx_in,
                  //               idx_out, lila::real(val));

                  fill(idx_out, idx_in, val);
                }
              }

              ++idx_dn;
            }
          }

          ++idx_up;
        }
      }

      // Apply hoppings on upspins
      if ((hop.type() == "HOP") || (hop.type() == "HOPUP")) {
        idx_t idx_up = 0;
        for (bit_t ups : block.reps_up_) {

	  if (popcnt(ups & flipmask) != 1) {
	    ++idx_up;
	    continue;
	  }

          idx_t up_offset_in = block.up_offsets_[idx_up];
          auto const &dnss_in = block.dns_for_up_rep(ups);
          auto const &norms_in = block.norms_for_up_rep(ups);

          bit_t ups_flip = ups ^ flipmask;
          idx_t idx_up_flip = block.index_up(ups_flip);
          idx_t up_offset_out = block.up_offsets_[idx_up_flip];
          bit_t ups_flip_rep = block.reps_up_[idx_up_flip];
          auto [sym_lower, sym_upper] = block.sym_limits_up(ups_flip);
          auto const &dnss_out = block.dns_for_up_rep(ups_flip_rep);
          auto const &norms_out = block.norms_for_up_rep(ups_flip_rep);

          double fermi_hop = popcnt(ups & spacemask) & 1 ? -1. : 1.;

          // trivial up-stabilizer (likely)
          if (sym_upper - sym_lower == 1) {
            int sym = block.syms_up_[sym_lower];
            bool fermi_up = block.fermi_bool_up(sym, ups_flip);
            idx_t idx_dn = 0;
            for (bit_t dns : dnss_in) {
              idx_t idx_in = up_offset_in + idx_dn;
              bit_t dns_rep = group_action.apply(sym, dns);
              idx_t idx_dn_flip = block.lintable_dns_.index(dns_rep);
              idx_t idx_out = up_offset_out + idx_dn_flip;
              double fermi = fermi_hop;
              if (fermi_up != block.fermi_bool_dn(sym, dns))
                fermi = -fermi;

              // Complex conjugate t if necessary
              if constexpr (is_complex<coeff_t>()) {
                val = -(gbit(dns, s1) ? t : tconj) * fermi *
                      irrep.character(sym) * norms_out[idx_dn_flip] /
                      norms_in[idx_dn];
              } else {
                val = -t * fermi * lila::real(irrep.character(sym)) *
                      norms_out[idx_dn_flip] / norms_in[idx_dn];
              }


                // idx_t idx_out2 = block.index(ups_flip, dns);
                // int n_sites = block.n_sites();
                // lila::Log.out("tri s1: {} s2: {} {};{} -> {};{} -> {};{} idx_in: "
                //               "{}, idx_out: {}, val: {:.4f}",
                //               s1, s2, bits_to_string(ups, n_sites),
                //               bits_to_string(dns, n_sites),
                //               bits_to_string(ups_flip, n_sites),
                //               bits_to_string(dns, n_sites),
                //               bits_to_string(ups_flip_rep, n_sites),
                //               bits_to_string(dns_rep, n_sites), idx_in, idx_out,
                //               lila::real(val));
                // assert(idx_out == idx_out2);

              fill(idx_out, idx_in, val);
              ++idx_dn;
            }
            // non-trivial up-stabilizer (unlikely)
          } else {
            idx_t idx_dn = 0;
            for (bit_t dns : dnss_in) {
              idx_t idx_in = up_offset_in + idx_dn;

              // Determine dns representative
              bit_t rep_dns = std::numeric_limits<bit_t>::max();
              int rep_sym = 0;
              for (idx_t sym_idx = sym_lower; sym_idx < sym_upper; ++sym_idx) {
                int sym = block.syms_up_[sym_idx];
                bit_t tdns = group_action.apply(sym, dns);
                if (tdns < rep_dns) {
                  rep_dns = tdns;
                  rep_sym = sym;
                }
              }
              bool fermi_up = block.fermi_bool_up(rep_sym, ups_flip);

              // Find index of rep dns
              auto it =
                  std::lower_bound(dnss_out.begin(), dnss_out.end(), rep_dns);

              if ((it != dnss_out.end()) && (*it == rep_dns)) {
                idx_t idx_dn_flip = std::distance(dnss_out.begin(), it);
                idx_t idx_out = up_offset_out + idx_dn_flip;
                double fermi = fermi_hop;
                if (fermi_up != block.fermi_bool_dn(rep_sym, dns))
                  fermi = -fermi;

                // Complex conjugate t if necessary
                if constexpr (is_complex<coeff_t>()) {
                  val = -(gbit(dns, s1) ? t : tconj) * fermi *
                        irrep.character(rep_sym) * norms_out[idx_dn_flip] /
                        norms_in[idx_dn];
                } else {
                  val = -t * fermi * lila::real(irrep.character(rep_sym)) *
                        norms_out[idx_dn_flip] / norms_in[idx_dn];
                }

                // // Debug
                // idx_t idx_out2 = block.index(ups_flip, dns);
                // int n_sites = block.n_sites();
                // lila::Log.out("non s1: {} s2: {} {};{} -> {};{} -> {};{} idx_in: "
                //               "{}, idx_out: {}/{}, val: {:.4f}",
                //               s1, s2, bits_to_string(ups, n_sites),
                //               bits_to_string(dns, n_sites),
                //               bits_to_string(ups_flip, n_sites),
                //               bits_to_string(dns, n_sites),
                //               bits_to_string(ups_flip_rep, n_sites),
                //               bits_to_string(rep_dns, n_sites), idx_in, idx_out, idx_out2,
                //               lila::real(val));
                // assert(idx_out == idx_out2);

                fill(idx_out, idx_in, val);
              }
              ++idx_dn;
            }
          }
          ++idx_up;
        }
      }
    }
  }
}
} // namespace hydra::electronterms
