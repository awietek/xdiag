#pragma once

#include <xdiag/blocks/electron.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis {

template <typename coeff_t>
void dispatch_matrix(OpSum const &ops, Electron const &block_in,
                     Electron const &block_out, coeff_t *mat, int64_t m);

} // namespace xdiag::basis
