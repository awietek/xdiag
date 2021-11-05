#pragma once

#include <lila/utils/logger.h>

<<<<<<< HEAD
#include <hydra/bitops/bitops.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/tj/tj_utils.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
=======
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/bitops/bitops.h>
>>>>>>> e87af3037790df1901285da9ed2ed0f3b66ef29d

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
      }
