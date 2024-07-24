#pragma once

#include <xdiag/blocks/tj_distributed.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis::tj_distributed {

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, tJDistributed const &block_in,
                    arma::Col<coeff_t> const &vec_in,
                    tJDistributed const &block_out,
                    arma::Col<coeff_t> &vec_out);

} // namespace xdiag::basis::tj_distributed
