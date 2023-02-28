#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/electron/electron.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void apply_ising_no_sym(Bond const &bond, Indexing &&indexing, Filler fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined() && (bond.type() == "ISING"));
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  coeff_t J = bond.coupling<coeff_t>();
  int s1 = std::min(bond.site(0), bond.site(1));
  int s2 = std::max(bond.site(0), bond.site(1));

  // Set values for same/diff (tJ block definition)
  coeff_t val_same = J / 4.;
  coeff_t val_diff = -J / 4.;

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

} // namespace hydra::electron
