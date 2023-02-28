#pragma once

namespace hydra::electron {

template <typename bit_t, typename coeff_t, class Indexing, class Fill>
void apply_number_sym(Bond const &bond, Indexing &&indexing, Fill fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined() &&
         ((bond.type() == "NUMBERUP") || (bond.type() == "NUMBERDN")));
  assert(bond.size() == 1);

  coeff_t mu = bond.coupling<coeff_t>();
  int s = bond.site(0);
  std::string type = bond.type();

  bit_t mask = (bit_t)1 << s;

  if (type == "NUMBERUP") {

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
  } else if (type == "NUMBERDN") {

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
  }
}

} // namespace hydra::electron
