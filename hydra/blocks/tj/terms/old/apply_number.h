#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void apply_number(Bond const &bond, Indexing &&indexing, Filler &&fill) {
  using namespace combinatorics;

  assert(bond.coupling_defined());
  assert(bond.type_defined() &&
         ((bond.type() == "NUMBERUP") || (bond.type() == "NUMBERDN")));
  assert(bond.size() == 1);

  int n_sites = indexing.n_sites();
  int n_up = indexing.n_up();
  int n_spins = indexing.n_spins();
  int n_holes = indexing.n_holes();

  coeff_t mu = bond.coupling<coeff_t>();
  int s = bond.site(0);
  std::string type = bond.type();

  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;
  bit_t mask = (bit_t)1 << s;
  idx_t idx = 0;
  for (bit_t holes : Combinations<bit_t>(n_sites, n_holes)) {

    // A hole is present at site -> continue
    if (holes & mask) {
      idx += binomial(n_spins, n_up);
      continue;
    }

    bit_t not_holes = (~holes) & sitesmask;
    if (type == "NUMBERUP") {
      for (bit_t spins : Combinations<bit_t>(n_spins, n_up)) {
        bit_t ups_dns = bitops::deposit(spins, not_holes);
        if (mask & ups_dns) {
          fill(idx, idx, mu);
        }
        ++idx;
      }
    } else if (type == "NUMBERDN") {
      for (bit_t spins : Combinations<bit_t>(n_spins, n_up)) {
        bit_t ups_dns = bitops::deposit(spins, not_holes);
        if (mask & (~ups_dns)) {
          fill(idx, idx, mu);
        }
        ++idx;
      }
    }
  }
}

} // namespace hydra::tj
