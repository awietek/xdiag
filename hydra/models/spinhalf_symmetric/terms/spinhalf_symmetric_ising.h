#pragma once
#include <tuple>

#include <lila/utils/logger.h>

#include <hydra/common.h>
#include <hydra/models/spinhalf_symmetric/spinhalf_symmetric.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <hydra/models/utils/model_utils.h>
#include <hydra/utils/bitops.h>

namespace hydra::terms::spinhalf_symmetric {

template <class bit_t, class Filler, class GroupAction>
void do_ising_symmetric(BondList const &bonds, Couplings const &couplings,
                        SpinhalfSymmetric<bit_t, GroupAction> const &block,
                        Filler &&fill) {

  auto ising = bonds.bonds_of_type("HEISENBERG") +
               bonds.bonds_of_type("ISING") + bonds.bonds_of_type("HB");

  for (auto bond : ising) {

    if (bond.size() != 2)
      lila::Log.err("Error computing SpinhalfSymmetric Ising: "
                    "bond must have exactly two sites defined");

    if (utils::coupling_is_non_zero(bond, couplings)) {

      std::string coupling = bond.coupling();
      double J = lila::real(couplings[coupling]);

      // Set values for same/diff
      std::string type = bond.type();
      double val_same = J / 4.;
      double val_diff = -J / 4.;

      int s1 = bond.site(0);
      int s2 = bond.site(1);
      if (s1 == s2)
        lila::Log.err("Error computing SpinhalfSymmetric Ising: "
                      "operator acting on twice the same site");
      bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

      for (idx_t idx = 0; idx < block.size(); ++idx) {
        bit_t state = block.indexing_.state(idx);
        if (utils::popcnt(state & mask) & 1) {
          fill(idx, idx, val_diff);
        } else {
          fill(idx, idx, val_same);
        }
      }
    }
  }
}

} // namespace hydra::terms::spinhalf_symmetric
