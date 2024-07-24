#pragma once

#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class Basis,
          class Fill>
void apply_ising(Op const &op, Basis &&basis, Fill &&fill) {
  assert(op.type() == "ISING");
  assert(op.size() == 2);
  assert(sites_disjoint(op));

  Coupling cpl = op.coupling();
  assert(cpl.isexplicit() && !cpl.ismatrix());
  coeff_t J = cpl.as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];

  // Set values for same/diff (tJ block definition)
  coeff_t val_same = J / 4.;
  coeff_t val_diff = -J / 4.;

  // bitmasks for fast evaluations
  bit_t s1mask = (bit_t)1 << s1;
  bit_t s2mask = (bit_t)1 << s2;
  bit_t mask = s1mask | s2mask;

  if constexpr (symmetric) {

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (int64_t idx_ups = 0; idx_ups < basis.n_rep_ups(); ++idx_ups) {
      bit_t ups = basis.rep_ups(idx_ups);
      auto dnss = basis.dns_for_ups_rep(ups);
      int64_t idx = basis.ups_offset(idx_ups);

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
      auto ups_and_idces = basis.states_indices_ups_thread();
#else
    auto ups_and_idces = basis.states_indices_ups();
#endif

      for (auto [up, idx_up] : ups_and_idces) {

        int64_t idx = idx_up * basis.size_dns();

        if ((up & mask) == mask) { // both spins pointing up
          for (bit_t dn : basis.states_dns()) {
            if (!(dn & mask)) {
              fill(idx, idx, val_same);
            }
            ++idx;
          }
        } else if (up & s1mask) { // s1 is pointing up
          for (bit_t dn : basis.states_dns()) {
            if ((dn & mask) == s2mask) {
              fill(idx, idx, val_diff);
            }
            ++idx;
          }
        } else if (up & s2mask) { // s2 is pointing up
          for (bit_t dn : basis.states_dns()) {
            if ((dn & mask) == s1mask) {
              fill(idx, idx, val_diff);
            }
            ++idx;
          }

        } else { // no upspins
          for (bit_t dn : basis.states_dns()) {
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

} // namespace xdiag::basis::electron
