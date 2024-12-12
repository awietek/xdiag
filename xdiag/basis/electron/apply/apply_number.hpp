#pragma once

namespace xdiag::basis::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class Basis,
          class Fill>
void apply_number(Coupling const &cpl, Op const &op, Basis &&basis,
                  Fill fill) try {
  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = (bit_t)1 << s;
  std::string type = op.type();

  if (type == "NUP") {

    if constexpr (symmetric) {
      int64_t idx = 0;
      for (int64_t idx_ups = 0; idx_ups < basis.n_rep_ups(); ++idx_ups) {
        bit_t ups = basis.rep_ups(idx_ups);
        auto dnss = basis.dns_for_ups_rep(ups);
        int64_t size_dns = dnss.size();
        if (ups & mask) { // check whether bit at s is set
          int64_t end = idx + size_dns;
          for (; idx < end; ++idx) {
            fill(idx, idx, mu);
          }
        } else {
          idx += size_dns;
        }
      }

    } else { // if not symmetric

#ifdef _OPENMP
#pragma omp parallel
      {
        auto ups_and_idces = basis.states_indices_ups_thread();
#else
      auto ups_and_idces = basis.states_indices_ups();
#endif
        for (auto [up, idx_up] : ups_and_idces) {
          int64_t idx = idx_up * basis.size_dns();

          int64_t size_dns = basis.size_dns();
          if (up & mask) { // check whether bit at s is set
            int64_t end = idx + size_dns;
            for (; idx < end; ++idx) {
              fill(idx, idx, mu);
            }
          } else {
            idx += size_dns;
          }
        }
#ifdef _OPENMP
      }
#endif
    }
  } else if (type == "NDN") {

    if constexpr (symmetric) {
      int64_t idx = 0;
      for (int64_t idx_ups = 0; idx_ups < basis.n_rep_ups(); ++idx_ups) {
        bit_t ups = basis.rep_ups(idx_ups);
        auto dnss = basis.dns_for_ups_rep(ups);
        for (bit_t dns : dnss) {
          if (dns & mask) {
            fill(idx, idx, mu);
          }
          ++idx;
        }
      }

    } else { // if not symmetric
#ifdef _OPENMP
#pragma omp parallel
      {
        auto ups_and_idces = basis.states_indices_ups_thread();
#else
      auto ups_and_idces = basis.states_indices_ups();
#endif

        for (auto [up, idx_up] : ups_and_idces) {
          (void)up;
          int64_t idx = idx_up * basis.size_dns();

          for (bit_t dn : basis.states_dns()) {
            if (dn & mask) {
              fill(idx, idx, mu);
            }
            ++idx;
          }
        }

#ifdef _OPENMP
      }
#endif
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::electron
