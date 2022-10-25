#pragma once

#include <hydra/blocks/tj/terms/apply_terms.h>
#include <hydra/common.h>
#include <hydra/indexing/indexing_variants.h>
#include <hydra/utils/logger.h>

namespace hydra::tj {

template <typename coeff_t, class Fill>
void apply_terms_dispatch(BondList const &bonds,
                          indexing::tJIndexing const &indexing_in,
                          indexing::tJIndexing const &indexing_out,
                          Fill &&fill) {
  using namespace indexing::tj;

  std::visit(
      overloaded{
          // uint16_t
          [&](IndexingNp<uint16_t> const &idx_in,
              IndexingNp<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSymmetricNp<uint16_t> const &idx_in,
              IndexingSymmetricNp<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          // uint32_t
          [&](IndexingNp<uint32_t> const &idx_in,
              IndexingNp<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSymmetricNp<uint32_t> const &idx_in,
              IndexingSymmetricNp<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          // uint64_t
          [&](IndexingNp<uint64_t> const &idx_in,
              IndexingNp<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSymmetricNp<uint64_t> const &idx_in,
              IndexingSymmetricNp<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](auto const &idx_in, auto const &idx_out) {
            Log.err("Error in apply_terms_dispatch: Invalid Indexing");
          }},
      indexing_in, indexing_out);
}

} // namespace hydra::tj
