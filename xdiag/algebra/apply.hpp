#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

void apply(OpSum const &ops, State const &v, State &w,
           double zero_precision = 1e-12);

template <typename coeff_t>
void apply(OpSum const &ops, Block const &block_in,
           arma::Col<coeff_t> const &vec_in, Block const &block_out,
           arma::Col<coeff_t> &vec_out, double zero_precision = 1e-12);

} // namespace xdiag
