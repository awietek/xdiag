#pragma once

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/blocks/blocks.h>
#include <xdiag/operators/bondlist.h>
#include <xdiag/states/state.h>

#include <xdiag/blocks/electron/electron_apply.h>
#include <xdiag/blocks/spinhalf/spinhalf_apply.h>
#include <xdiag/blocks/tj/tj_apply.h>
#include <xdiag/blocks/tj_distributed/tj_distributed_apply.h>

namespace xdiag {

void apply(BondList const &bonds, State const &v, State &w);
  
template <typename coeff_t>
void apply(BondList const &bonds, block_variant_t const &block_in,
           arma::Col<coeff_t> const &vec_in, block_variant_t const &block_out,
           arma::Col<coeff_t> &vec_out);


} // namespace xdiag
