#pragma once

#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/blocks/spinhalf/spinhalf.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/bondlist.hpp>

namespace xdiag {

template <typename coeff_t>
void apply(BondList const &bonds, Spinhalf const &block_in,
           arma::Col<coeff_t> const &vec_in, Spinhalf const &block_out,
           arma::Col<coeff_t> &vec_out);

} // namespace xdiag
