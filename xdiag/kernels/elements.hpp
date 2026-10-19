// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <functional>

#include <xdiag/kernels/compile.hpp>

namespace xdiag::kernels {

// A type-erased callback receiving one matrix element (row, col, value).
template <typename coeff_t>
using fill_t = std::function<void(int64_t, int64_t, coeff_t)>;

// Internal (not part of the public API).
//
// Calls fill(row, col, val) once for every nonzero matrix element
// val = ops(row, col) of the compiled operator (row in block_out, col in
// block_in). Rows are split over the OpenMP threads: all calls for one row come
// from one thread, in a fixed order, so fill may write to anything indexed by
// row without synchronisation. Throws for a distributed block.
template <typename coeff_t>
void elements(CompiledOpSum const &ops, fill_t<coeff_t> const &fill);

} // namespace xdiag::kernels
