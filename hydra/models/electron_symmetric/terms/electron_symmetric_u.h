#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/models/electron_symmetric/electron_symmetric.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>
#include <hydra/symmetries/symmetry_utils.h>

namespace hydra::electron {

template <class bit_t, class SymmetryGroup, class Filler>
void do_U_symmetric(Couplings const &couplings,
                    ElectronSymmetric<bit_t, SymmetryGroup> const &block,
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

} // namespace hydra::electron
