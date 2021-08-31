#pragma once

#include <lila/utils/logger.h>

#include <hydra/common.h>
#include <hydra/models/spinhalf_symmetric/spinhalf_symmetric.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/utils/model_utils.h>

namespace hydra::spinhalfterms {

template <class bit_t, class coeff_t, class Filler, class GroupAction>
void do_exchange_symmetric(BondList const &bonds, Couplings const &couplings,
                           SpinhalfSymmetric<bit_t, GroupAction> const &block,
                           Filler &&fill) {

  auto exchange = bonds.bonds_of_type("HEISENBERG") +
                  bonds.bonds_of_type("EXCHANGE") + bonds.bonds_of_type("HB");
  auto const &group_action = block.group_action();
  auto const &irrep = block.irrep();

  for (auto bond : exchange) {

    if (bond.size() != 2)
      lila::Log.err("Error computing Spinhalf Exchange: "
                    "bond must have exactly two sites defined");

    if (utils::coupling_is_non_zero(bond, couplings)) {
      std::string coupling = bond.coupling();
      double J = lila::real(couplings[coupling]); // TODO: could be complex

      // Set values for same/diff (tJ model definition)
      double Jhalf = J / 2.;
      int s1 = bond.site(0);
      int s2 = bond.site(1);
      if (s1 == s2) {
        // TODO: this case could actually be implemented
        lila::Log.err("Error computing Spinhalf Exchange: "
                      "operator acting on twice the same site");
      }

      bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
      idx_t idx_in = 0;
      for (bit_t state : block.states_) {

        if (utils::popcnt(state & mask) & 1) {

          bit_t new_state = state ^ mask;
          auto [rep, rep_sym] = group_action.representative_index(new_state);
          idx_t idx_out = block.index(rep);

	  // if new state has non-zero norm, compute element and fill
          if (idx_out != invalid_index) {
	    coeff_t val = Jhalf * complex_to<coeff_t>(irrep.character(rep_sym)) *
	      block.norms_[idx_out] / block.norms_[idx_in];
	    fill(idx_out, idx_in, val);
	  }

        }
        ++idx_in;
      }
    }
  }
}

} // namespace hydra::spinhalfterms
