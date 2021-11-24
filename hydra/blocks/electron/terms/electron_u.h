#pragma once

#include <lila/utils/logger.h>

#include <hydra/common.h>

#include <hydra/blocks/utils/block_utils.h>

#include <hydra/bitops/bitops.h>

#include <hydra/indexing/electron/electron_indexing.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <hydra/combinatorics/combinations.h>

namespace hydra::terms {

template <class bit_t, class coeff_t, class Filler>
void electron_U(Couplings const &couplings,
                indexing::ElectronIndexing<bit_t> const &indexing,
                Filler &&fill) {
  using bitops::popcnt;
  using combinatorics::Combinations;

  int n_sites = indexing.n_sites();
  int n_up = indexing.n_up();
  int n_dn = indexing.n_dn();

  if (couplings.defined("U")) {

    coeff_t U = utils::get_coupling<coeff_t>(couplings, "U");

    if (!lila::close(U, 0.)) {
      idx_t idx = 0;
      for (bit_t up : Combinations<bit_t>(n_sites, n_up)) {
        for (bit_t dn : Combinations<bit_t>(n_sites, n_dn)) {
          coeff_t val = U * (double)popcnt(up & dn);
          fill(idx, idx, val);
          ++idx;
        }
      }
    }
  }
}

} // namespace hydra::terms
