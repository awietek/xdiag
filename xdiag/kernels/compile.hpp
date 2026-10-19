// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::kernels {

// Internal (not part of the public API).
//
// An OpSum prepared for the kernels, mapping block_in onto block_out. The
// kernels enumerate their input states in parallel, each thread owning a
// contiguous range, and write only at the input index (the gather fill). To
// apply ops we therefore hand them hc(ops), whose input states are the rows of
// ops: every output element is written by one thread in a fixed order, without
// atomics, so the result is bitwise identical for any number of threads.
struct CompiledOpSum {
  // hc(ops), normal-ordered with the block's implementation algebra. For a
  // distributed block ops itself: its MPI kernels normal-order on their own.
  OpSum ops;
  Block block_in;  // domain of the original ops
  Block block_out; // codomain of the original ops
};

// Checks that ops maps block_in onto block_out and compiles it. Iterative
// algorithms compile once and apply the result on every iteration.
CompiledOpSum compile(OpSum const &ops, Block const &block_in,
                      Block const &block_out);

} // namespace xdiag::kernels
