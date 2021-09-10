#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/models/electron_symmetric/electron_symmetric.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/symmetry_utils.h>
#include <hydra/utils/bitops.h>

namespace hydra::electronterms {

template <class bit_t, class GroupAction, class Filler>
void do_U_symmetric(Couplings const &couplings,
                    ElectronSymmetric<bit_t, GroupAction> const &block,
                    Filler &&fill) {
  using utils::popcnt;
  if (couplings.defined("U")) {

    if (!couplings.is_real("U")) {
      lila::Log.err("Error computing ElectronSymmetric matrix/apply: "
                    "Hubbard U must be a real number");
    }

    double U = couplings.real("U");
    if (!lila::close(U, 0.)) {

      idx_t idx = 0;
      for (bit_t ups : block.reps_up_) {
	auto const &dnss = block.dns_for_up_rep(ups);
	for (bit_t dns : dnss) {
          double val = U * (double)popcnt(ups & dns);
          fill(idx, idx, val);
	  ++idx;
        }
      }
    }
  }
}

} // namespace hydra::electronterms
