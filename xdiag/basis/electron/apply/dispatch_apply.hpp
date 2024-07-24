#pragma once

#include <xdiag/blocks/electron.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis::electron {

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, Electron const &block_in,
                    arma::Col<coeff_t> const &vec_in, Electron const &block_out,
                    arma::Col<coeff_t> &vec_out);

} // namespace xdiag::basis::electron
