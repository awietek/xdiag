#pragma once

#include <lila/utils/logger.h>

#include <hydra/common.h>
#include <hydra/blocks/electron_symmetric_simple/electron_symmetric_simple.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/symmetry_operations.h>
#include <hydra/bitops/bitops.h>

namespace hydra::terms::electron_symmetric_simple {

template <class bit_t, class coeff_t, class GroupAction, class Filler>
void do_hopping_symmetric(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetricSimple<bit_t, GroupAction> const &block, Filler &&fill) {
  using bitops::gbit;
  using bitops::popcnt;

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
    if (couplings.defined(cpl) && !lila::close(couplings[cpl], (complex)0.)) {
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

        for (auto [up, lower_upper] : block.ups_lower_upper_) {
          idx_t lower = lower_upper.first;
          idx_t upper = lower_upper.second;
          auto begin_up = block.dns_.begin() + lower;
          auto end_up = block.dns_.begin() + upper;

          std::vector<int> stable_syms = group_action.stabilizer_symmetries(up);

          // trivial stabilizer of ups -> no representative lookup -> faster
          if (stable_syms.size() == 1) {
            for (idx_t idx = lower; idx < upper; ++idx) {
              bit_t dn = block.dn(idx);

              // If hopping is possible ...
              if (popcnt(dn & flipmask) == 1) {
                bit_t dn_flip = dn ^ flipmask;

                // Look for index os dn flip
                auto it = std::lower_bound(begin_up, end_up, dn_flip);

                // if a state has been found
                if ((it != end_up) && (*it == dn_flip)) {
                  idx_t idx_out = std::distance(begin_up, it) + lower;
                  double fermi_hop = popcnt(dn & spacemask) & 1 ? -1. : 1.;

                  // Complex conjugate t if necessary
                  if constexpr (is_complex<coeff_t>()) {
                    val = -(gbit(dn, s1) ? t : tconj) * fermi_hop *
                          block.norm(idx_out) / block.norm(idx);
                  } else {
                    val =
                        -t * fermi_hop * block.norm(idx_out) / block.norm(idx);
                  }
                  fill(idx_out, idx, val);
                }
              }
            }
          } else { // non-trivial stabilizer of ups
            for (idx_t idx = lower; idx < upper; ++idx) {
              bit_t dn = block.dn(idx);

              // If hopping is possible ...
              if (popcnt(dn & flipmask) == 1) {
                bit_t dn_flip = dn ^ flipmask;

                // Determine the dn representative from stable symmetries
                bit_t dn_flip_rep = dn_flip;
                int dn_flip_rep_sym = stable_syms[0];
                for (auto const &stable_sym : stable_syms) {
                  bit_t trans = group_action.apply(stable_sym, dn_flip);
                  if (trans < dn_flip_rep) {
                    dn_flip_rep = trans;
                    dn_flip_rep_sym = stable_sym;
                  }
                }

                // Look for index os dn flip
                auto it = std::lower_bound(begin_up, end_up, dn_flip_rep);

                // if a state has been found
                if ((it != end_up) && (*it == dn_flip_rep)) {
                  idx_t idx_out = std::distance(begin_up, it) + lower;

                  double fermi_hop = popcnt(dn & spacemask) & 1 ? -1. : 1.;
                  double fermi_up =
                      group_action.fermi_sign(dn_flip_rep_sym, up);
                  double fermi_dn =
                      group_action.fermi_sign(dn_flip_rep_sym, dn_flip);

                  // Complex conjugate t if necessary
                  if constexpr (is_complex<coeff_t>()) {
                    val = -(gbit(dn, s1) ? t : tconj) * fermi_hop * fermi_up *
                          fermi_dn * irrep.character(dn_flip_rep_sym) *
                          block.norm(idx_out) / block.norm(idx);
                  } else {
                    val = -t * fermi_hop * fermi_up * fermi_dn *
                          lila::real(irrep.character(dn_flip_rep_sym)) *
                          block.norm(idx_out) / block.norm(idx);
                  }

                  fill(idx_out, idx, val);
                }
              }
            } // non-trivial stabilizer of ups

          } // for (auto [up, lower_upper]
        }
      } //     if ((hop.type() == "HOP") || (hop.type() == "HOPDN")) {

      // Apply hoppings on upspins
      if ((hop.type() == "HOP") || (hop.type() == "HOPUP")) {

        for (auto [dn, lower_upper] : block.dns_lower_upper_) {
          idx_t lower = lower_upper.first;
          idx_t upper = lower_upper.second;
          auto begin_dn = block.ups_.begin() + lower;
          auto end_dn = block.ups_.begin() + upper;
          std::vector<int> stable_syms = group_action.stabilizer_symmetries(dn);

          // trivial stabilizer of dns -> no representative lookup -> faster
          if (stable_syms.size() == 1) {
            for (idx_t idx = lower; idx < upper; ++idx) {
              bit_t up = block.up(idx);

              // If hopping is possible ...
              if (popcnt(up & flipmask) == 1) {
                bit_t up_flip = up ^ flipmask;
                idx_t idx_in = block.index_switch_to_index(idx);
                coeff_t chi_switch_in =
                    complex_to<coeff_t>(block.character_switch_[idx]);

                // Look for index of up_flip
                auto it = std::lower_bound(begin_dn, end_dn, up_flip);

                // if a state has been found
                if ((it != end_dn) && (*it == up_flip)) {
                  idx_t idx_out_switch = std::distance(begin_dn, it) + lower;
                  idx_t idx_out = block.index_switch_to_index(idx_out_switch);

                  double fermi_hop = popcnt(up & spacemask) & 1 ? -1. : 1.;
                  coeff_t chi_switch_out = complex_to<coeff_t>(
                      block.character_switch_[idx_out_switch]);

                  // Complex conjugate t if necessary
                  if constexpr (is_complex<coeff_t>()) {
                    val = -(gbit(up, s1) ? t : tconj) * fermi_hop *
                          lila::conj(chi_switch_in) * chi_switch_out *
                          block.norm(idx_out) / block.norm(idx_in);
                  } else {
                    val = -t * fermi_hop * chi_switch_in * chi_switch_out *
                          block.norm(idx_out) / block.norm(idx_in);
                  }
                  fill(idx_out, idx_in, val);
                }
              }
            }
          } else { // non-trivial stabilizer of ups
            for (idx_t idx = lower; idx < upper; ++idx) {
              bit_t up = block.up(idx);

              // If hopping is possible ...
              if (popcnt(up & flipmask) == 1) {
                bit_t up_flip = up ^ flipmask;
                idx_t idx_in = block.index_switch_to_index(idx);
                coeff_t chi_switch_in =
                    complex_to<coeff_t>(block.character_switch_[idx]);

                // Determine the up representative from stable symmetries
                bit_t up_flip_rep = up_flip;
                int up_flip_rep_sym = stable_syms[0];
                for (auto const &stable_sym : stable_syms) {
                  bit_t trans = group_action.apply(stable_sym, up_flip);
                  if (trans < up_flip_rep) {
                    up_flip_rep = trans;
                    up_flip_rep_sym = stable_sym;
                  }
                }

                // Look for index os dn flip
                auto it = std::lower_bound(begin_dn, end_dn, up_flip_rep);

                // if a state has been found
                if ((it != end_dn) && (*it == up_flip_rep)) {
                  idx_t idx_out_switch = std::distance(begin_dn, it) + lower;
                  idx_t idx_out = block.index_switch_to_index(idx_out_switch);

                  double fermi_hop = popcnt(up & spacemask) & 1 ? -1. : 1.;
                  double fermi_up =
                      group_action.fermi_sign(up_flip_rep_sym, up_flip);
                  double fermi_dn =
                      group_action.fermi_sign(up_flip_rep_sym, dn);
                  coeff_t chi_switch_out = complex_to<coeff_t>(
                      block.character_switch_[idx_out_switch]);

                  // Complex conjugate t if necessary
                  if constexpr (is_complex<coeff_t>()) {
                    val = -(gbit(up, s1) ? t : tconj) * fermi_hop * fermi_up *
                          fermi_dn * irrep.character(up_flip_rep_sym) *
                          lila::conj(chi_switch_in) * chi_switch_out *
                          block.norm(idx_out) / block.norm(idx_in);
                  } else {
                    val = -t * fermi_hop * fermi_up * fermi_dn *
                          lila::real(irrep.character(up_flip_rep_sym)) *
                          chi_switch_in * chi_switch_out * block.norm(idx_out) /
                          block.norm(idx_in);
                  }
                  fill(idx_out, idx_in, val);
                }
              }
            }
          } // non-trivial stabilizer of ups
        }
      }
    }
  }
}
} // namespace hydra::terms::electron_symmetric_simple
