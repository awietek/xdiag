#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

#include <xdiag/blocks/electron/electron_apply.hpp>
#include <xdiag/blocks/spinhalf/spinhalf_apply.hpp>
#include <xdiag/blocks/tj/tj_apply.hpp>
#include <xdiag/blocks/tj_distributed/tj_distributed_apply.hpp>

namespace xdiag {

void apply(OpSum const &ops, State const &v, State &w);

template <typename coeff_t>
void apply(OpSum const &ops, block_variant_t const &block_in,
           arma::Col<coeff_t> const &vec_in, block_variant_t const &block_out,
           arma::Col<coeff_t> &vec_out);

} // namespace xdiag
