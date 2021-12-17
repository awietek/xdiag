#pragma once

#include <lila/utils/logger.h>

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/electron/electron.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/operators/operator_utils.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void electron_ising(BondList const &bonds, Couplings const &couplings,
                    Indexing &&indexing, Filler &&fill) {
  using bitops::gbit;
  using bitops::popcnt;

  auto clean_bonds = utils::clean_bondlist(
      bonds, couplings,
      {"HEISENBERG", "HB", "ISING", "TJHEISENBERG", "TJISING", "TJHB"}, 2);

  for (auto bond : clean_bonds) {

    std::string type = bond.type();
    std::string cpl = bond.coupling();

    utils::check_sites_disjoint(bond);
    int s1 = std::min(bond.site(0), bond.site(1));
    int s2 = std::max(bond.site(0), bond.site(1));

    coeff_t J = utils::get_coupling<coeff_t>(couplings, cpl);

    // Set values for same/diff (tJ block definition)
    coeff_t val_same, val_diff;
    if ((type == "HEISENBERG") || (type == "ISING") || (type == "HB")) {
      val_same = J / 4.;
      val_diff = -J / 4.;
    } else if ((type == "TJHEISENBERG") || (type == "TJISING") ||
               (type == "TJHB")) {
      val_same = 0.;
      val_diff = -J / 2.;
    }

    bit_t s1mask = (bit_t)1 << s1;
    bit_t s2mask = (bit_t)1 << s2;
    bit_t mask = s1mask | s2mask;

    idx_t idx = 0;
    for (auto up : indexing.states_ups()) {

      if ((up & mask) == mask) { // both spins pointing up
        for (bit_t dn : indexing.states_dns()) {
          if (!(dn & mask)) {
            fill(idx, idx, val_same);
          }
          ++idx;
        }
      } else if (up & s1mask) { // s1 is pointing up
        for (bit_t dn : indexing.states_dns()) {
          if ((dn & mask) == s2mask) {
            fill(idx, idx, val_diff);
          }
          ++idx;
        }
      } else if (up & s2mask) { // s2 is pointing up
        for (bit_t dn : indexing.states_dns()) {
          if ((dn & mask) == s1mask) {
            fill(idx, idx, val_diff);
          }
          ++idx;
        }

      } else { // no upspins
        for (bit_t dn : indexing.states_dns()) {
          if ((dn & mask) == mask) {
            fill(idx, idx, val_same);
          }
          ++idx;
        }
      }

    } // for (auto up : ...)
  }
}

} // namespace hydra::terms
