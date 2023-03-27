#pragma once

namespace hydra::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class Indexing,
          class Fill>
void apply_number(Bond const &bond, Indexing &&indexing, Fill fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert(bond.size() == 1);

  std::string type = bond.type();
  assert((type == "NUMBERUP") || (type == "NUMBERDN"));

  coeff_t mu = bond.coupling<coeff_t>();
  int s = bond.site(0);

  bit_t mask = (bit_t)1 << s;

  if (type == "NUMBERUP") {

    if constexpr (symmetric) {
      idx_t idx = 0;
      for (idx_t idx_ups = 0; idx_ups < indexing.n_rep_ups(); ++idx_ups) {
        bit_t ups = indexing.rep_ups(idx_ups);
        auto dnss = indexing.dns_for_ups_rep(ups);
        idx_t size_dns = dnss.size();
        if (ups & mask) { // check whether bit at s is set
          idx_t end = idx + size_dns;
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
        auto ups_and_idces = indexing.states_indices_ups_thread();
#else
      auto ups_and_idces = indexing.states_indices_ups();
#endif
        for (auto [up, idx_up] : ups_and_idces) {
          idx_t idx = idx_up * indexing.size_dns();

          idx_t size_dns = indexing.size_dns();
          if (up & mask) { // check whether bit at s is set
            idx_t end = idx + size_dns;
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
  } else if (type == "NUMBERDN") {

    if constexpr (symmetric) {
      idx_t idx = 0;
      for (idx_t idx_ups = 0; idx_ups < indexing.n_rep_ups(); ++idx_ups) {
        bit_t ups = indexing.rep_ups(idx_ups);
        auto dnss = indexing.dns_for_ups_rep(ups);
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
        auto ups_and_idces = indexing.states_indices_ups_thread();
#else
      auto ups_and_idces = indexing.states_indices_ups();
#endif

        for (auto [up, idx_up] : ups_and_idces) {

          idx_t idx = idx_up * indexing.size_dns();

          for (bit_t dn : indexing.states_dns()) {
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
}

} // namespace hydra::electron
