#pragma once

#include <lila/utils/logger.h>

#include <hydra/common.h>
#include <hydra/bitops/bitops.h>

#include <hydra/blocks/blocks.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/symmetries/symmetry_utils.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms::electron_symmetric {

template <class bit_t, class coeff_t, class GroupAction, class Filler>
void do_exchange_symmetric(BondList const &bonds, Couplings const &couplings,
                           ElectronSymmetric<bit_t, GroupAction> const &block,
                           Filler &&fill) {
  using bitops::gbit;
  using bitops::popcnt;

  auto const &group_action = block.group_action();
  auto const &irrep = block.irrep();

  auto exchange = bonds.bonds_of_type("HEISENBERG") +
                  bonds.bonds_of_type("EXCHANGE") + bonds.bonds_of_type("HB");
  auto exchange_tj =
      bonds.bonds_of_type("TJHEISENBERG") + bonds.bonds_of_type("TJHB");

  for (auto bond : exchange + exchange_tj) {

    if (bond.size() != 2)
      lila::Log.err("Error computing Electron exchange: "
                    "bonds must have exactly two sites defined");

    if (!utils::coupling_is_zero(bond, couplings)) {
      int s1 = std::min(bond.site(0), bond.site(1));
      int s2 = std::max(bond.site(0), bond.site(1));
      bit_t s1mask = (bit_t)1 << s1;
      bit_t s2mask = (bit_t)1 << s2;
      bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
      bit_t spacemask = ((1 << (s2 - s1 - 1)) - 1)
                        << (s1 + 1); // used to determine Fermi sign

      std::string cpl = bond.coupling();
      double J = lila::real(couplings[cpl]);
      double Jhalf = J / 2.;
      coeff_t val;

      // Loop over all up configurations
      idx_t idx_up = 0;
      for (bit_t ups : block.reps_up_) {

        // Exchange gives zero continue
        if ((popcnt(ups & mask) == 2) || (popcnt(ups & mask) == 0)) {
          ++idx_up;
          continue;
        }

        // Get the flip masks for up and down spins
        bit_t dns_mask = ((ups & mask) == s1mask) ? s2mask : s1mask;

        // Get the origin dns for this ups
        idx_t up_offset_in = block.up_offsets_[idx_up];
        auto const &dnss_in = block.dns_for_up_rep(ups);
        auto const &norms_in = block.norms_for_up_rep(ups);

        // Get flipped ups, its representative, symmetries, and target dns
        bit_t ups_flip = ups ^ mask;
        idx_t idx_up_flip = block.index_up(ups_flip);
        idx_t up_offset_out = block.up_offsets_[idx_up_flip];
        bit_t ups_flip_rep = block.reps_up_[idx_up_flip];
        auto [sym_lower, sym_upper] = block.sym_limits_up(ups_flip);
        assert(sym_upper - sym_lower > 0);
        auto const &dnss_out = block.dns_for_up_rep(ups_flip_rep);
        auto const &norms_out = block.norms_for_up_rep(ups_flip_rep);

        // flag whether or not we have a sign change
        bool fermi_up = !(bool)(popcnt(ups & spacemask) & 1);

        // trivial up-stabilizer (likely)
        if (sym_upper - sym_lower == 1) {
          int sym = block.syms_up_[sym_lower];
          fermi_up ^= block.fermi_bool_up(sym, ups_flip);

          idx_t idx_dn = 0;
          for (bit_t dns : dnss_in) {

            // If  dns can be raised
            if ((dns & mask) == dns_mask) {
              idx_t idx_in = up_offset_in + idx_dn;
              bit_t dns_flip = dns ^ mask;
              bit_t dns_rep = group_action.apply(sym, dns_flip);
              idx_t idx_dn_flip = block.lintable_dns_.index(dns_rep);
              idx_t idx_out = up_offset_out + idx_dn_flip;
              bool fermi = fermi_up;
              fermi ^= (bool)(popcnt(dns & spacemask) & 1);
              fermi ^= block.fermi_bool_dn(sym, dns_flip);

              if constexpr (is_complex<coeff_t>()) {
                val = Jhalf * irrep.character(sym) * norms_out[idx_dn_flip] /
                      norms_in[idx_dn];
              } else {
                val = Jhalf * lila::real(irrep.character(sym)) *
                      norms_out[idx_dn_flip] / norms_in[idx_dn];
              }

              // idx_t idx_out2 = block.index(ups_flip, dns_flip);
              // int n_sites = block.n_sites();
              // lila::Log.out("tri s1: {} s2: {} {};{} -> {};{} -> {};{}
              // idx_in: "
              //               "{}, idx_out: {}/{}, val: {:.4f}",
              //               s1, s2, bits_to_string(ups, n_sites),
              //               bits_to_string(dns, n_sites),
              //               bits_to_string(ups_flip, n_sites),
              //               bits_to_string(dns_flip, n_sites),
              //               bits_to_string(ups_flip_rep, n_sites),
              //               bits_to_string(dns_rep, n_sites), idx_in,
              //               idx_out, idx_out2, lila::real(val));
              // assert(idx_out == idx_out2);

              // lila::Log.out("{} {} {} {} {}", (bool)(popcnt(ups & spacemask)
              // & 1),
              //               block.fermi_bool_up(sym, ups_flip),
              //               (bool)(popcnt(dns & spacemask) & 1),
              //               block.fermi_bool_dn(sym, dns_flip), fermi);

              if (fermi)
                fill(idx_out, idx_in, -val);
              else
                fill(idx_out, idx_in, val);
            }
            ++idx_dn;
          }
          // non-trivial up-stabilizer (unlikely)
        } else {
          idx_t idx_dn = 0;
          for (bit_t dns : dnss_in) {

            // If  dns can be raised
            if ((dns & mask) == dns_mask) {

              idx_t idx_in = up_offset_in + idx_dn;
              bit_t dns_flip = dns ^ mask;

              // Determine dns representative
              bit_t rep_dns = std::numeric_limits<bit_t>::max();
              int rep_sym = 0;
              for (idx_t sym_idx = sym_lower; sym_idx < sym_upper; ++sym_idx) {
                int sym = block.syms_up_[sym_idx];
                bit_t tdns = group_action.apply(sym, dns_flip);
                if (tdns < rep_dns) {
                  rep_dns = tdns;
                  rep_sym = sym;
                }
              }

              // Find index of rep dns
              auto it =
                  std::lower_bound(dnss_out.begin(), dnss_out.end(), rep_dns);

              if ((it != dnss_out.end()) && (*it == rep_dns)) {
                idx_t idx_dn_flip = std::distance(dnss_out.begin(), it);
                idx_t idx_out = up_offset_out + idx_dn_flip;
                bool fermi = fermi_up ^ block.fermi_bool_up(rep_sym, ups_flip);
                fermi ^= (bool)(popcnt(dns & spacemask) & 1);
                fermi ^= block.fermi_bool_dn(rep_sym, dns_flip);

                if constexpr (is_complex<coeff_t>()) {
                  val = Jhalf * irrep.character(rep_sym) *
                        norms_out[idx_dn_flip] / norms_in[idx_dn];
                } else {
                  val = Jhalf * lila::real(irrep.character(rep_sym)) *
                        norms_out[idx_dn_flip] / norms_in[idx_dn];
                }

                // // Debug
                // idx_t idx_out2 = block.index(ups_flip, dns_flip);
                // idx_t idx_in2 = block.index(ups, dns);

                // int n_sites = block.n_sites();
                // lila::Log.out(
                //     "non s1: {} s2: {} {};{} -> {};{} -> {};{} idx_in: "
                //     "{}, idx_out: {}/{}, val: {:.4f}",
                //     s1, s2, bits_to_string(ups, n_sites),
                //     bits_to_string(dns, n_sites),
                //     bits_to_string(ups_flip, n_sites),
                //     bits_to_string(dns_flip, n_sites),
                //     bits_to_string(ups_flip_rep, n_sites),
                //     bits_to_string(rep_dns, n_sites), idx_in, idx_out,
                //     idx_out2, lila::real(val));
                // assert(idx_in == idx_in2);
                // assert(idx_out == idx_out2);

                if (fermi)
                  fill(idx_out, idx_in, -val);
                else
                  fill(idx_out, idx_in, val);
              }
            }
            ++idx_dn;
          }
        }
        ++idx_up;
      }
    }
  }
}

} // namespace hydra::terms::electron_symmetric
