#pragma once

#include <hydra/common.h>

#include <hydra/blocks/utils/block_utils.h>

#include <hydra/bitops/bitops.h>

#include <hydra/indexing/electron/electron_symmetric_indexing.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>


namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Filler>
void electron_symmetric_ising(
    BondList const &bonds, Couplings const &couplings,
    indexing::ElectronSymmetricIndexing<bit_t> const &indexing, Filler &&fill) {
  auto clean_bonds = utils::clean_bondlist(
      bonds, couplings,
      {"HEISENBERG", "HB", "ISING", "TJHEISENBERG", "TJISING", "TJHB"}, 2);

  for (auto bond : clean_bonds) {

    std::string type = bond.type();
    std::string cpl = bond.coupling();

    utils::check_sites_disjoint(bond);
    int s1 = bond[0];
    int s2 = bond[1];

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

    // bitmasks for fast evaluations
    bit_t s1mask = (bit_t)1 << s1;
    bit_t s2mask = (bit_t)1 << s2;
    bit_t mask = s1mask | s2mask;

    idx_t idx = 0;
    for (idx_t idx_ups = 0; idx_ups < indexing.n_rep_ups(); ++idx_ups) {
      bit_t ups = indexing.rep_ups(idx_ups);
      auto dnss = indexing.dns_for_ups_rep(ups);

      if ((ups & mask) == mask) { // both spins pointing up
        for (bit_t dns : dnss) {
          if (!(dns & mask))
            fill(idx, idx, val_same);
          ++idx;
        }
      } else if (ups & s1mask) { // s1 is pointing up
        for (bit_t dns : dnss) {
          if ((dns & mask) == s2mask)
            fill(idx, idx, val_diff);
          ++idx;
        }
      } else if (ups & s2mask) { // s2 is pointing up
        for (bit_t dns : dnss) {
          if ((dns & mask) == s1mask)
            fill(idx, idx, val_diff);
          ++idx;
        }
      } else { // no upspins
        for (bit_t dns : dnss) {
          if ((dns & mask) == mask)
            fill(idx, idx, val_same);
          ++idx;
        }
      }
    }
  }
}

} // namespace hydra::terms::electron_symmetric
