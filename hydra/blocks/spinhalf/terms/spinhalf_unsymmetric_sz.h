#pragma once
#include <hydra/common.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/operators/operator_utils.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void spinhalf_unsymmetric_sz(BondList const &bonds, Couplings const &couplings,
                             Indexing &&indexing, Filler &&fill) {
  auto clean_bonds = utils::clean_bondlist(bonds, couplings, {"SZ"}, 1);

  for (auto bond : clean_bonds) {

    utils::check_sites_disjoint(bond);
    int s = bond[0];
    bit_t mask = ((bit_t)1 << s);

    std::string cpl = bond.coupling();
    coeff_t J = utils::get_coupling<coeff_t>(couplings, cpl);
    coeff_t val_up = J / 2.;
    coeff_t val_dn = -J / 2.;

    idx_t idx = 0;
    for (auto spins : indexing.states()) {
      if (spins & mask)
        fill(idx, idx, val_up);
      else
        fill(idx, idx, val_dn);

      ++idx;
    }
  }
}

} // namespace hydra::terms
