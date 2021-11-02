#pragma once

#include <lila/utils/logger.h>

#include <hydra/common.h>
#include <hydra/blocks/spinhalf_symmetric/spinhalf_symmetric.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra::terms::spinhalf_symmetric {

template <class bit_t, class coeff_t, class Filler, class GroupAction>
void do_exchange_symmetric(BondList const &bonds, Couplings const &couplings,
                           SpinhalfSymmetric<bit_t, GroupAction> const &block,
                           Filler &&fill) {

  auto exchange = bonds.bonds_of_type("HEISENBERG") +
                  bonds.bonds_of_type("EXCHANGE") + bonds.bonds_of_type("HB");

  for (auto bond : exchange) {

    if (bond.size() != 2)
      lila::Log.err("Error computing Spinhalf Exchange: "
                    "bond must have exactly two sites defined");

    if (utils::coupling_is_non_zero(bond, couplings)) {

      // Get couplings and define flip mask
      std::string coupling = bond.coupling();
      double J = lila::real(couplings[coupling]); // TODO: could be complex
      double Jhalf = J / 2.;
      int s1 = bond.site(0);
      int s2 = bond.site(1);
      if (s1 == s2) {
        lila::Log.err("Error computing Spinhalf Exchange: "
                      "operator acting on twice the same site");
      }
      bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

      // Run through all states and apply bond
      for (idx_t idx_in = 0; idx_in < block.size(); ++idx_in) {
        bit_t state_in = block.indexing_.state(idx_in);

        if (bitops::popcnt(state_in & mask) & 1) { // if spins are flippable
          bit_t state_out = state_in ^ mask;
          idx_t idx_out = block.indexing_.index(state_out);

          // if new state has non-zero norm, compute element and fill
          if (idx_out != invalid_index) {
            auto norm_out =
                complex_to<coeff_t>(block.indexing_.norm(state_out));
            auto norm_in = complex_to<coeff_t>(block.indexing_.norm(state_in));
            fill(idx_out, idx_in, Jhalf * norm_out / norm_in);
          }
        }
      } // for (idx_t idx_in; ...)
    }
  }
}

} // namespace hydra::terms::spinhalf_symmetric
