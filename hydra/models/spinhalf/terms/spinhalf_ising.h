#pragma once
#include <tuple>

#include <lila/utils/logger.h>

#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/models/electron/electron.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>

namespace hydra::spinhalfterms {
  
template <class bit_t, class Filler>
void do_ising(BondList const &bonds, Couplings const &couplings,
	      Spinhalf<bit_t> const &block, Filler &&fill) {

  auto ising = bonds.bonds_of_type("HEISENBERG") +
               bonds.bonds_of_type("ISING") + bonds.bonds_of_type("HB");

  for (auto bond : ising) {

    if (bond.size() != 2)
      lila::Log.err("Error computing Spinhalf Ising: "
                   "bond must have exactly two sites defined");

    std::string coupling = bond.coupling();
    if (couplings.defined(coupling) &&
        !lila::close(couplings[coupling], (complex)0.)) {

      double J = lila::real(couplings[coupling]);

      // Set values for same/diff (tJ model definition)
      std::string type = bond.type();
      double val_same = J / 4.;
      double val_diff = -J / 4.;

      int s1 = bond.site(0);
      int s2 = bond.site(1);
      if (s1 == s2)
	lila::Log.err("Error computing Spinhalf Ising: "
		     "operator acting on twice the same site");
      bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

      int n_sites = block.n_sites();
      int n_up = block.n_up();
      idx_t idx = 0;
      for (auto spins : Combinations<bit_t>(n_sites, n_up)) {

        if (utils::popcnt(spins & mask) & 1)
          fill(idx, idx, val_diff);
        else
          fill(idx, idx, val_same);

        ++idx;
      }
    }
  }
}

} // namespace hydra::spinhalfterms
