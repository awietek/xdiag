#pragma once

#include <xdiag/operators/opsum.hpp>
#include <xdiag/blocks/spinhalf.hpp>

namespace xdiag::basis::spinhalf {

template <typename coeff_t>
void dispatch_matrix(OpSum const &ops, Spinhalf const &block_in,
                     Spinhalf const &block_out, coeff_t *mat, int64_t m);

} // namespace xdiag::basis::spinhalf
