#pragma once

#include <hydra/common.h>
#include <hydra/models/electron.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>

namespace hydra::electron {

template <class bit_t, class Filler>
void do_U(Couplings const &couplings, Electron<bit_t> const &block,
          Filler &&filler) {
  using utils::popcnt;

  if (couplings.defined("U")) {

    int n_sites = block.n_sites();
    int n_up = block.n_up();
    int n_dn = block.n_dn();

    if (!couplings.is_real("U")) {
      HydraLog.err("Error creating Electron matrix: "
                   "Hubbard U must be a real number");
    }

    double U = couplings.real("U");
    if (!lila::close(U, 0.)) {
      idx_t idx = 0;
      for (auto up : Combinations(n_sites, n_up)) {
        for (auto dn : Combinations(n_sites, n_dn)) {
          filler(idx, U, up, dn);
          ++idx;
        }
      }
    }
  }
}

} // namespace hydra::electron
