#pragma once
#ifdef XDIAG_USE_MPI

#include <xdiag/blocks/tj_distributed/tj_distributed.hpp>
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

template <typename coeff_t>
void apply(OpSum const &ops, tJDistributed const &block_in,
           arma::Col<coeff_t> const &vec_in, tJDistributed const &block_out,
           arma::Col<coeff_t> &vec_out);

} // namespace xdiag
#endif
