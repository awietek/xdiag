#pragma once
#include <string>

#include <hydra/common.h>

#include <hydra/combinatorics/combinations.h>

#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/utils/block_utils.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/operators/operator_utils.h>

#include <hydra/bitops/bitops.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class IndexingIn, class IndexingOut,
          class Filler>
void spinhalf_unsymmetric_spsm(BondList const &bonds,
                               Couplings const &couplings,
                               IndexingIn &&indexing_in,
                               IndexingOut &&indexing_out, Filler &&fill,
                               std::string spsm) {
  assert((spsm == "S+") || (spsm == "S-"));
  auto clean_bonds = utils::clean_bondlist(bonds, couplings, {spsm}, 1);

  for (auto bond : clean_bonds) {
    
    int s = bond[0];
    bit_t mask = ((bit_t)1 << s);
    std::string cpl = bond.coupling();
    coeff_t J = utils::get_coupling<coeff_t>(couplings, cpl);

    idx_t idx_in = 0;

    if (spsm == "S+") {
      for (auto spins_in : indexing_in.states()) {
        if (!(spins_in & mask)) {
          bit_t spins_out = spins_in | mask;
          idx_t idx_out = indexing_out.index(spins_out);
          fill(idx_out, idx_in, J);
        }
        ++idx_in;
      }
    } else if (spsm == "S-") {
      for (auto spins_in : indexing_in.states()) {
        if (spins_in & mask) {
          bit_t spins_out = spins_in ^ mask;
          idx_t idx_out = indexing_out.index(spins_out);
          fill(idx_out, idx_in, J);
        }
        ++idx_in;
      }
    }

  } // for (auto bond : clean_bonds)

}

} // namespace hydra::terms
