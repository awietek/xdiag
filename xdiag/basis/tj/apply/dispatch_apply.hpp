#pragma once

#include <xdiag/blocks/tj.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis::tj {

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, tJ const &block_in,
                    arma::Col<coeff_t> const &vec_in, tJ const &block_out,
                    arma::Col<coeff_t> &vec_out);

} // namespace xdiag::basis::tj
