#pragma once

#include <lila/utils/logger.h>

#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/blocks/electron/electron.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>

namespace hydra::terms::electron {

template <class bit_t, class Filler>
void do_U(Couplings const &couplings, Electron<bit_t> const &block,
          Filler &&fill) {
  using combinatorics::Combinations;
  using bitops::popcnt;

  if (couplings.defined("U")) {

    int n_sites = block.n_sites();
    int n_up = block.n_up();
    int n_dn = block.n_dn();

    if (!couplings.is_real("U")) {
      lila::Log.err("Error creating Electron matrix: "
                    "Hubbard U must be a real number");
    }

    double U = couplings.real("U");
    if (!lila::close(U, 0.)) {
      idx_t idx = 0;
      for (bit_t up : Combinations<bit_t>(n_sites, n_up)) {
        for (bit_t dn : Combinations<bit_t>(n_sites, n_dn)) {
          double val = U * popcnt(up & dn);
          fill(idx, idx, val);
          ++idx;
        }
      }
    }
  }
}

} // namespace hydra::terms::electron
