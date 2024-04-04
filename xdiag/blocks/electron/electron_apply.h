#pragma once

#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/common.h>
#include <xdiag/operators/bondlist.h>

#include <xdiag/blocks/electron/electron.h>


namespace xdiag {

template <typename coeff_t>
void apply(BondList const &bonds, Electron const &block_in,
           arma::Col<coeff_t> const &vec_in, Electron const &block_out,
           arma::Col<coeff_t> &vec_out);

} // namespace xdiag
