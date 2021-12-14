#pragma once

#include <hydra/common.h>

#include <hydra/blocks/utils/block_utils.h>

#include <hydra/bitops/bitops.h>

#include <hydra/indexing/electron/electron_symmetric_indexing.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void electron_symmetric_U(Couplings const &couplings, Indexing &&indexing,
                          Filler &&fill) {

  if (couplings.defined("U")) {

    coeff_t U = utils::get_coupling<coeff_t>(couplings, "U");

    if (!lila::close(U, 0.)) {

      idx_t idx = 0;
      for (idx_t idx_ups = 0; idx_ups < indexing.n_rep_ups(); ++idx_ups) {
        bit_t ups = indexing.rep_ups(idx_ups);
        auto dnss = indexing.dns_for_ups_rep(ups);
        for (bit_t dns : dnss) {
          coeff_t val = U * (double)bitops::popcnt(ups & dns);
          fill(idx, idx, val);
          ++idx;
        }
      }
    }
  }
}

} // namespace hydra::terms
