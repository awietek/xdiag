#pragma once

#include <variant>

#include <hydra/common.h>
#include <hydra/utils/logger.h>

#include <hydra/indexing/indexing_variants.h>

#include <hydra/blocks/spinhalf/terms/apply_terms.h>

namespace hydra::terms::spinhalf {

using namespace hydra::indexing::spinhalf;

template <typename bit_t, typename coeff_t, class Fill>
void apply_terms_dispatch(BondList const &bonds, Couplings const &couplings,
                          Indexing<bit_t> const &indexing_in,
                          Indexing<bit_t> const &indexing_out, Fill &&fill) {

  std::visit(
      overloaded{[&](IndexingSz<bit_t> const &idx_in,
                     IndexingSz<bit_t> const &idx_out) {
                   terms::spinhalf::apply_terms<bit_t, coeff_t, false>(
                       bonds, couplings, idx_in, idx_out, fill);
                 },
                 [&](IndexingNoSz<bit_t> const &idx_in,
                     IndexingNoSz<bit_t> const &idx_out) {
                   terms::spinhalf::apply_terms<bit_t, coeff_t, false>(
                       bonds, couplings, idx_in, idx_out, fill);
                 },
                 [&](IndexingSymmetricSz<bit_t> const &idx_in,
                     IndexingSymmetricSz<bit_t> const &idx_out) {
                   apply_terms<bit_t, coeff_t, true>(bonds, couplings, idx_in,
                                                     idx_out, fill);
                 },
                 [&](IndexingSymmetricNoSz<bit_t> const &idx_in,
                     IndexingSymmetricNoSz<bit_t> const &idx_out) {
                   apply_terms<bit_t, coeff_t, true>(bonds, couplings, idx_in,
                                                     idx_out, fill);
                 },
                 [&](IndexingSublattice<bit_t, 1> const &idx_in,
                     IndexingSublattice<bit_t, 1> const &idx_out) {
                   apply_terms<bit_t, coeff_t, true>(bonds, couplings, idx_in,
                                                     idx_out, fill);
                 },
                 [&](IndexingSublattice<bit_t, 2> const &idx_in,
                     IndexingSublattice<bit_t, 2> const &idx_out) {
                   apply_terms<bit_t, coeff_t, true>(bonds, couplings, idx_in,
                                                     idx_out, fill);
                 },
                 [&](IndexingSublattice<bit_t, 3> const &idx_in,
                     IndexingSublattice<bit_t, 3> const &idx_out) {
                   apply_terms<bit_t, coeff_t, true>(bonds, couplings, idx_in,
                                                     idx_out, fill);
                 },
                 [&](IndexingSublattice<bit_t, 4> const &idx_in,
                     IndexingSublattice<bit_t, 4> const &idx_out) {
                   apply_terms<bit_t, coeff_t, true>(bonds, couplings, idx_in,
                                                     idx_out, fill);
                 },
                 [&](IndexingSublattice<bit_t, 5> const &idx_in,
                     IndexingSublattice<bit_t, 5> const &idx_out) {
                   apply_terms<bit_t, coeff_t, true>(bonds, couplings, idx_in,
                                                     idx_out, fill);
                 },
                 [&](auto const &idx_in, auto const &idx_out) {
                   Log.err("Error in apply_terms_dispatch: Invalid Indexing");
                 }},
      indexing_in, indexing_out);
}

} // namespace hydra::terms::spinhalf
