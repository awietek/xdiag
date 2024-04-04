#pragma once

#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/common.h>
#include <xdiag/operators/bondlist.h>

#include <xdiag/blocks/tj/tj.h>

namespace xdiag {

template <typename coeff_t>
void apply(BondList const &bonds, tJ const &block_in,
           arma::Col<coeff_t> const &vec_in, tJ const &block_out,
           arma::Col<coeff_t> &vec_out);

} // namespace xdiag
