#pragma once

#include <xdiag/blocks/electron/electron.hpp>
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

template <typename coeff_t>
void apply(OpSum const &ops, Electron const &block_in,
           arma::Col<coeff_t> const &vec_in, Electron const &block_out,
           arma::Col<coeff_t> &vec_out, double zero_precision = 1e-12);

} // namespace xdiag
