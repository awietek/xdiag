#pragma once

#include <xdiag/blocks/tj/tj.hpp>
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/bondlist.hpp>

namespace xdiag {

template <typename coeff_t>
void apply(BondList const &bonds, tJ const &block_in,
           arma::Col<coeff_t> const &vec_in, tJ const &block_out,
           arma::Col<coeff_t> &vec_out, double zero_precision = 1e-12);

} // namespace xdiag
