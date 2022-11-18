#pragma once

#include <hydra/common.h>
#include <hydra/utils/logger.h>

#include <hydra/indexing/indexing_variants.h>

#include <hydra/blocks/spinhalf/terms/apply_terms.h>

namespace hydra::spinhalf {

template <typename coeff_t, class Fill>
void apply_terms_dispatch(BondList const &bonds,
                          indexing::SpinhalfIndexing const &indexing_in,
                          indexing::SpinhalfIndexing const &indexing_out,
                          Fill &&fill) {
  using namespace indexing::spinhalf;

  std::visit(
      variant::overloaded{
          // uint16_t
          [&](IndexingSz<uint16_t> const &idx_in,
              IndexingSz<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingNoSz<uint16_t> const &idx_in,
              IndexingNoSz<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSymmetricSz<uint16_t> const &idx_in,
              IndexingSymmetricSz<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSymmetricNoSz<uint16_t> const &idx_in,
              IndexingSymmetricNoSz<uint16_t> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint16_t, 1> const &idx_in,
              IndexingSublattice<uint16_t, 1> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint16_t, 2> const &idx_in,
              IndexingSublattice<uint16_t, 2> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint16_t, 3> const &idx_in,
              IndexingSublattice<uint16_t, 3> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint16_t, 4> const &idx_in,
              IndexingSublattice<uint16_t, 4> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint16_t, 5> const &idx_in,
              IndexingSublattice<uint16_t, 5> const &idx_out) {
            apply_terms<uint16_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },

          // uint32_t
          [&](IndexingSz<uint32_t> const &idx_in,
              IndexingSz<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingNoSz<uint32_t> const &idx_in,
              IndexingNoSz<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSymmetricSz<uint32_t> const &idx_in,
              IndexingSymmetricSz<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSymmetricNoSz<uint32_t> const &idx_in,
              IndexingSymmetricNoSz<uint32_t> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint32_t, 1> const &idx_in,
              IndexingSublattice<uint32_t, 1> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint32_t, 2> const &idx_in,
              IndexingSublattice<uint32_t, 2> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint32_t, 3> const &idx_in,
              IndexingSublattice<uint32_t, 3> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint32_t, 4> const &idx_in,
              IndexingSublattice<uint32_t, 4> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint32_t, 5> const &idx_in,
              IndexingSublattice<uint32_t, 5> const &idx_out) {
            apply_terms<uint32_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },

          // uint64_t
          [&](IndexingSz<uint64_t> const &idx_in,
              IndexingSz<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingNoSz<uint64_t> const &idx_in,
              IndexingNoSz<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, false>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSymmetricSz<uint64_t> const &idx_in,
              IndexingSymmetricSz<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSymmetricNoSz<uint64_t> const &idx_in,
              IndexingSymmetricNoSz<uint64_t> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint64_t, 1> const &idx_in,
              IndexingSublattice<uint64_t, 1> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint64_t, 2> const &idx_in,
              IndexingSublattice<uint64_t, 2> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint64_t, 3> const &idx_in,
              IndexingSublattice<uint64_t, 3> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint64_t, 4> const &idx_in,
              IndexingSublattice<uint64_t, 4> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },
          [&](IndexingSublattice<uint64_t, 5> const &idx_in,
              IndexingSublattice<uint64_t, 5> const &idx_out) {
            apply_terms<uint64_t, coeff_t, true>(bonds, idx_in, idx_out, fill);
          },

          [&](auto const &idx_in, auto const &idx_out) {
            Log.err("Error in apply_terms_dispatch: Invalid Indexing");
          }},
      indexing_in, indexing_out);
}

} // namespace hydra::spinhalf
