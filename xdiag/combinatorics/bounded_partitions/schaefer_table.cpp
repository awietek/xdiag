// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "schaefer_table.hpp"

#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/pack_unpack.hpp>
#include <xdiag/combinatorics/bounded_partitions/count_bounded_partitions.hpp>
#include <xdiag/math/ipow.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::combinatorics {

template <typename bitarray_t>
SchaeferTable<bitarray_t>::SchaeferTable(int64_t n, int64_t total,
                                         int64_t d) try
    : bounded_partitions_(n, total, d), n_fast_((n + 1) / 2),
      n_slow_(n / 2), d_(d),
      fast_rank_table_(math::ipow(d, (n + 1) / 2), 0),
      slow_offset_table_(math::ipow(d, n / 2), 0) {

  // ------------------------------------------------------------------
  // Build fast_rank_table_:
  //   For each k_fast in [0, n_fast_*(d-1)], iterate
  //   BoundedPartitions(n_fast_, k_fast, d) in rlex order and store
  //   the local rank 0,1,2,... in fast_rank_table_[fast_pack].
  //   k_fast is uniquely determined by fast_pack (it is the digit sum
  //   of fast_pack in base d), so there are no conflicts.
  // ------------------------------------------------------------------
  for (int64_t k_fast = 0; k_fast <= n_fast_ * (d - 1); ++k_fast) {
    int64_t local_rank = 0;
    for (auto pattern : BoundedPartitions<bitarray_t>(n_fast_, k_fast, d)) {
      int64_t fp = bits::pack(pattern, d, n_fast_);
      fast_rank_table_[fp] = local_rank++;
    }
  }

  // ------------------------------------------------------------------
  // Build slow_offset_table_:
  //   Iterate slow_pack = 0, 1, ..., d^n_slow_ - 1.  Ascending
  //   slow_pack corresponds to rlex order on the slow subsystem (since
  //   the little-endian base-d encoding has the lowest-index slow
  //   slot as the least significant digit = fastest-varying in rlex).
  //   For each slow_pack, decode the digit sum (slow_sum), compute
  //   k_fast = total - slow_sum, and accumulate
  //   count_BP(n_fast_, k_fast, d) as the offset increment.
  // ------------------------------------------------------------------
  int64_t slow_table_size = math::ipow(d, n_slow_);
  int64_t offset = 0;
  for (int64_t sp = 0; sp < slow_table_size; ++sp) {
    slow_offset_table_[sp] = offset;
    // Decode digit sum of sp in base d
    int64_t tmp = sp, slow_sum = 0;
    for (int64_t i = 0; i < n_slow_; ++i) {
      slow_sum += tmp % d;
      tmp /= d;
    }
    int64_t k_fast = total - slow_sum;
    int64_t cnt = (k_fast >= 0 && k_fast <= n_fast_ * (d - 1))
                      ? count_bounded_partitions(n_fast_, k_fast, d)
                      : 0;
    offset += cnt;
  }
}
XDIAG_CATCH

template <typename bitarray_t> int64_t SchaeferTable<bitarray_t>::n() const {
  return bounded_partitions_.n();
}

template <typename bitarray_t>
int64_t SchaeferTable<bitarray_t>::total() const {
  return bounded_partitions_.total();
}

template <typename bitarray_t>
int64_t SchaeferTable<bitarray_t>::d() const {
  return bounded_partitions_.d();
}

template <typename bitarray_t> int64_t SchaeferTable<bitarray_t>::size() const {
  return bounded_partitions_.size();
}

template <typename bitarray_t>
int64_t SchaeferTable<bitarray_t>::bitwidth() const {
  return bounded_partitions_.bitwidth();
}

template <typename bitarray_t>
bitarray_t SchaeferTable<bitarray_t>::operator[](int64_t idx) const {
  return bounded_partitions_[idx];
}

template <typename bitarray_t>
BoundedPartitionsIterator<bitarray_t> SchaeferTable<bitarray_t>::begin() const {
  return bounded_partitions_.begin();
}

template <typename bitarray_t>
BoundedPartitionsIterator<bitarray_t> SchaeferTable<bitarray_t>::end() const {
  return bounded_partitions_.end();
}

template <typename bitarray_t>
bool SchaeferTable<bitarray_t>::operator==(
    SchaeferTable<bitarray_t> const &rhs) const {
  return bounded_partitions_ == rhs.bounded_partitions_;
}

template <typename bitarray_t>
bool SchaeferTable<bitarray_t>::operator!=(
    SchaeferTable<bitarray_t> const &rhs) const {
  return !operator==(rhs);
}

} // namespace xdiag::combinatorics

using namespace xdiag::bits;

#define INSTANTIATE_ST(BIT_T, NBITS)                                           \
  template class xdiag::combinatorics::SchaeferTable<                          \
      xdiag::bits::BitArray<BIT_T, NBITS>>;

#define INSTANTIATE_ST_ALL(BIT_T)                                              \
  INSTANTIATE_ST(BIT_T, 1)                                                     \
  INSTANTIATE_ST(BIT_T, 2)                                                     \
  INSTANTIATE_ST(BIT_T, 3)                                                     \
  INSTANTIATE_ST(BIT_T, 4)                                                     \
  INSTANTIATE_ST(BIT_T, 8)

// BEGIN_INSTANTIATION_GROUP(native)
INSTANTIATE_ST_ALL(uint16_t)
INSTANTIATE_ST_ALL(uint32_t)
INSTANTIATE_ST_ALL(uint64_t)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(bitset)
INSTANTIATE_ST_ALL(BitsetDynamic)
INSTANTIATE_ST_ALL(BitsetStatic1)
INSTANTIATE_ST_ALL(BitsetStatic2)
INSTANTIATE_ST_ALL(BitsetStatic4)
INSTANTIATE_ST_ALL(BitsetStatic8)
// END_INSTANTIATION_GROUP

#undef INSTANTIATE_ST_ALL
#undef INSTANTIATE_ST
