#pragma once

#include <xdiag/blocks/spinhalf_distributed.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis {

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, SpinhalfDistributed const &block_in,
                    arma::Col<coeff_t> const &vec_in,
                    SpinhalfDistributed const &block_out,
                    arma::Col<coeff_t> &vec_out);

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, SpinhalfDistributed const &block_in,
                    arma::Mat<coeff_t> const &vec_in,
                    SpinhalfDistributed const &block_out,
                    arma::Mat<coeff_t> &vec_out);
  
} // namespace xdiag::basis
