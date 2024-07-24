#pragma once

#include <xdiag/operators/opsum.hpp>
#include <xdiag/blocks/electron.hpp>

namespace xdiag::basis::electron {

template <typename coeff_t>
void dispatch_matrix(OpSum const &ops, Electron const &block_in,
                     Electron const &block_out, coeff_t *mat, int64_t m);

} // namespace xdiag::basis::electron
