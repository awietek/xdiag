#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class Indexing,
          class Fill>
void apply_ising(Bond const &bond, Indexing &&indexing, Fill &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined() && (bond.type() == "ISING"));
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  coeff_t J = bond.coupling<coeff_t>();
  int s1 = bond[0];
  int s2 = bond[1];

  // Set values for same/diff (tJ block definition)
  coeff_t val_same = J / 4.;
  coeff_t val_diff = -J / 4.;

  // bitmasks for fast evaluations
  bit_t s1mask = (bit_t)1 << s1;
  bit_t s2mask = (bit_t)1 << s2;
  bit_t mask = s1mask | s2mask;

  if constexpr (symmetric) {

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

  } else { // if not (symmetric)

#ifdef _OPENMP
#pragma omp parallel
    {
      auto ups_and_idces = indexing.states_indices_ups_thread();
#else
    auto ups_and_idces = indexing.states_indices_ups();
#endif

      for (auto [up, idx_up] : ups_and_idces) {

        idx_t idx = idx_up * indexing.size_dns();

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

#ifdef _OPENMP
    }
#endif
  }
}

} // namespace hydra::electron
