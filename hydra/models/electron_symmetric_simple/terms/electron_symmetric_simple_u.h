#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/models/electron_symmetric_simple/electron_symmetric_simple.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/symmetry_utils.h>
#include <hydra/utils/bitops.h>

namespace hydra::terms::electron_symmetric_simple {

template <class bit_t, class GroupAction, class Filler>
void do_U_symmetric(Couplings const &couplings,
                    ElectronSymmetricSimple<bit_t, GroupAction> const &block,
                    Filler &&fill) {
  using utils::popcnt;
  if (couplings.defined("U")) {

    if (!couplings.is_real("U")) {
      lila::Log.err("Error creating Electron matrix: "
                    "Hubbard U must be a real number");
    }

    double U = couplings.real("U");
    if (!lila::close(U, 0.)) {
      for (auto [up, lower_upper] : block.ups_lower_upper_) {
        idx_t lower = lower_upper.first;
        idx_t upper = lower_upper.second;
        for (idx_t idx = lower; idx < upper; ++idx) {
          bit_t dn = block.dn(idx);
          double val = U * popcnt(up & dn);
          fill(idx, idx, val);
          // mat(idx, idx) += val;
        }
      }
    }
  }
}

} // namespace hydra::terms::electron_symmetric_simple
