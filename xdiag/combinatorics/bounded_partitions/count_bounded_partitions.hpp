// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

namespace xdiag::combinatorics {

// ---------------------------------------------------------------------------
// Count weak compositions of total into n parts each in {0, ..., d-1}.
// Inclusion-exclusion: sum_{k=0}^{K} (-1)^k C(n,k) C(total-k*d+n-1, n-1)
// ---------------------------------------------------------------------------
int64_t count_bounded_partitions(int64_t n, int64_t total, int64_t d);

// ---------------------------------------------------------------------------
// Compute the sequence at linear index idx in rlex order.
// Processes slots from n-1 (most significant) down to 0: for each slot,
// enumerate values v = 0,1,... and subtract block counts until idx falls
// within the block for v.
// ---------------------------------------------------------------------------
template <typename bitarray_t>
bitarray_t nth_bounded_partition(int64_t n, int64_t total, int64_t d,
                                 int64_t idx);

// ---------------------------------------------------------------------------
// Compute the linear index (rank) of a sequence in rlex order.
// Inverse of nth_bounded_partition.
// ---------------------------------------------------------------------------
template <typename bitarray_t>
int64_t rank_bounded_partition(int64_t n, int64_t total, int64_t d,
                               bitarray_t seq);

} // namespace xdiag::combinatorics
