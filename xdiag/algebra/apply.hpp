#pragma once

#include <xdiag/common.hpp>

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/blocks/spinhalf_distributed.hpp>
#include <xdiag/blocks/tj_distributed.hpp>
#include <xdiag/blocks/electron_distributed.hpp>
#endif

namespace xdiag {

XDIAG_API State apply(Op const &op, State const &v);
XDIAG_API State apply(OpSum const &ops, State const &v);
XDIAG_API void apply(Op const &op, State const &v, State &w);
XDIAG_API void apply(OpSum const &ops, State const &v, State &w);

template <typename mat_t>
void apply(OpSum const &ops, Block const &block_in, mat_t const &mat_in,
           Block const &block_out, mat_t &mat_out);

template <typename mat_t, typename block_t>
void apply(OpSum const &ops, block_t const &block_in, mat_t const &mat_in,
           block_t const &block_out, mat_t &mat_out);

} // namespace xdiag
