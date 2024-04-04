#pragma once
#ifdef XDIAG_USE_MPI

#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/blocks/tj_distributed/tj_distributed.h>
#include <xdiag/common.h>
#include <xdiag/operators/bondlist.h>

namespace xdiag {

template <typename coeff_t>
void apply(BondList const &bonds, tJDistributed const &block_in,
           arma::Col<coeff_t> const &vec_in, tJDistributed const &block_out,
           arma::Col<coeff_t> &vec_out);

} // namespace xdiag
#endif
