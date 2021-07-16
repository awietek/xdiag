#pragma once

#include <lila/utils/logger.h>

#include <hydra/common.h>
#include <hydra/models/electron_symmetric/electron_symmetric.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>
#include <hydra/symmetries/symmetry_utils.h>

namespace hydra::electron {

template <class bit_t, class coeff_t, class SymmetryGroup, class Filler>
void do_hopping_symmetric(BondList const &bonds, Couplings const &couplings,
                          ElectronSymmetric<bit_t, SymmetryGroup> const &block,
                          Filler &&fill) {
  using utils::gbit;
  using utils::popcnt;

  auto symmetry_group = block.symmetry_group();
  auto irrep = block.irrep();

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

      // Apply hoppings on dnspins
      if ((hop.type() == "HOP") || (hop.type() == "HOPDN")) {

        for (auto [up, lower_upper] : block.ups_lower_upper_) {
          idx_t lower = lower_upper.first;
          idx_t upper = lower_upper.second;
          auto begin_up = block.dns_.begin() + lower;
          auto end_up = block.dns_.begin() + upper;

          std::vector<int> stable_syms =
	    utils::stabilizer_symmetries(up, symmetry_group);
          auto stable_group = symmetry_group.subgroup(stable_syms);
          auto stable_irrep = irrep.subgroup(stable_syms);

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

                  coeff_t t = (gbit(dn, s1)) ? couplings[cpl]
                                             : lila::conj(couplings[cpl]);
                  double fermi_hop = popcnt(dn & spacemask) & 1 ? -1. : 1.;
                  coeff_t val =
                      -t * fermi_hop * block.norm(idx_out) / block.norm(idx);

                  // mat(idx_out, idx) += val;
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
                auto [dn_flip_rep, dn_flip_rep_sym] =
                    stable_group.representative_index(dn_flip);

                // Look for index os dn flip
                auto it = std::lower_bound(begin_up, end_up, dn_flip_rep);

                // if a state has been found
                if ((it != end_up) && (*it == dn_flip_rep)) {
                  idx_t idx_out = std::distance(begin_up, it) + lower;

                  coeff_t t = (gbit(dn, s1)) ? couplings[cpl]
                                             : lila::conj(couplings[cpl]);
                  double fermi_hop = popcnt(dn & spacemask) & 1 ? -1. : 1.;
                  double fermi_up =
                      stable_group.fermi_sign(dn_flip_rep_sym, up);
                  double fermi_dn =
                      stable_group.fermi_sign(dn_flip_rep_sym, dn_flip);
                  coeff_t val = -t * fermi_hop * fermi_up * fermi_dn *
                                stable_irrep.character(dn_flip_rep_sym) *
                                block.norm(idx_out) / block.norm(idx);

                  // mat(idx_out, idx) += val;
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

          std::vector<int> stable_syms =
              utils::stabilizer_symmetries(dn, symmetry_group);
          auto stable_group = symmetry_group.subgroup(stable_syms);
          auto stable_irrep = irrep.subgroup(stable_syms);

          // trivial stabilizer of dns -> no representative lookup -> faster
          if (stable_syms.size() == 1) {
            for (idx_t idx = lower; idx < upper; ++idx) {
              bit_t up = block.up(idx);

              // If hopping is possible ...
              if (popcnt(up & flipmask) == 1) {
                bit_t up_flip = up ^ flipmask;
                idx_t idx_in = block.index_switch_to_index(idx);
                coeff_t chi_switch_in = block.character_switch_[idx];

                // Look for index of up_flip
                auto it = std::lower_bound(begin_dn, end_dn, up_flip);

                // if a state has been found
                if ((it != end_dn) && (*it == up_flip)) {
                  idx_t idx_out_switch = std::distance(begin_dn, it) + lower;
                  idx_t idx_out = block.index_switch_to_index(idx_out_switch);

                  coeff_t t = (gbit(up, s1)) ? couplings[cpl]
                                             : lila::conj(couplings[cpl]);
                  double fermi_hop = popcnt(up & spacemask) & 1 ? -1. : 1.;
                  coeff_t chi_switch_out =
                      block.character_switch_[idx_out_switch];
                  coeff_t val = -t * fermi_hop * lila::conj(chi_switch_in) *
                                chi_switch_out * block.norm(idx_out) /
                                block.norm(idx_in);
                  // mat(idx_out, idx_in) += val;
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
                coeff_t chi_switch_in = block.character_switch_[idx];

                // Determine the dn representative from stable symmetries
                auto [up_flip_rep, up_flip_rep_sym] =
                    stable_group.representative_index(up_flip);

                // Look for index os dn flip
                auto it = std::lower_bound(begin_dn, end_dn, up_flip_rep);

                // if a state has been found
                if ((it != end_dn) && (*it == up_flip_rep)) {
                  idx_t idx_out_switch = std::distance(begin_dn, it) + lower;
                  idx_t idx_out = block.index_switch_to_index(idx_out_switch);

                  coeff_t t = (gbit(up, s1)) ? couplings[cpl]
                                             : lila::conj(couplings[cpl]);
                  double fermi_hop = popcnt(up & spacemask) & 1 ? -1. : 1.;
                  double fermi_up =
                      stable_group.fermi_sign(up_flip_rep_sym, up_flip);
                  double fermi_dn =
                      stable_group.fermi_sign(up_flip_rep_sym, dn);
                  coeff_t chi_switch_out =
                      block.character_switch_[idx_out_switch];
                  coeff_t val = -t * fermi_hop * fermi_up * fermi_dn *
                                stable_irrep.character(up_flip_rep_sym) *
                                lila::conj(chi_switch_in) * chi_switch_out *
                                block.norm(idx_out) / block.norm(idx_in);

                  // mat(idx_out, idx_in) += val;
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
} // namespace hydra::electron
