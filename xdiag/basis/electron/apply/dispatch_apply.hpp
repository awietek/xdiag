#pragma once

#include <xdiag/blocks/electron.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis {

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, Electron const &block_in,
                    arma::Col<coeff_t> const &vec_in, Electron const &block_out,
                    arma::Col<coeff_t> &vec_out);

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, Electron const &block_in,
                    arma::Mat<coeff_t> const &mat_in, Electron const &block_out,
                    arma::Mat<coeff_t> &mat_out);

} // namespace xdiag::basis
