#pragma once

namespace hydra::electron {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void apply_number(Bond const &bond, Indexing &&indexing, Filler fill) {
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
    for (bit_t up : indexing.states_ups()) {

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
  } else if (type == "NUMBERDN") {

    idx_t idx = 0;
    idx_t size_ups = indexing.size_ups();
    for (idx_t up_idx = 0; up_idx < size_ups; ++up_idx) {
      for (bit_t dn : indexing.states_dns()) {
        if (dn & mask) {
          fill(idx, idx, mu);
        }
        ++idx;
      }
    }
  }
}

} // namespace hydra::electron
