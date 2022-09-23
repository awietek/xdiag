#pragma once

#include <hydra/blocks/tj/terms/apply_terms.h>
#include <hydra/common.h>
#include <hydra/indexing/indexing_variants.h>
#include <hydra/utils/logger.h>

namespace hydra::tj {

using namespace hydra::indexing::tj;

template <typename bit_t, typename coeff_t, class Fill>
void apply_terms_dispatch(BondList const &bonds,
                          Indexing<bit_t> const &indexing_in,
                          Indexing<bit_t> const &indexing_out, Fill &&fill) {

  std::visit(
      overloaded{
          [&](IndexingNp<bit_t> const &idx_in,
              IndexingNp<bit_t> const &idx_out) {
            apply_terms<bit_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSymmetricNp<bit_t> const &idx_in,
              IndexingSymmetricNp<bit_t> const &idx_out) {
            apply_terms<bit_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](auto const &idx_in, auto const &idx_out) {
            Log.err("Error in apply_terms_dispatch: Invalid Indexing");
          }},
      indexing_in, indexing_out);
}

} // namespace hydra::tj
