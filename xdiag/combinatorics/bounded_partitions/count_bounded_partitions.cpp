// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "count_bounded_partitions.hpp"

#include <algorithm>
#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/math/binomial.hpp>

namespace xdiag::combinatorics {

int64_t count_bounded_partitions(int64_t n, int64_t total, int64_t d) {
  if (n == 0)
    return (total == 0) ? 1 : 0;
  if (total < 0 || total > n * (d - 1))
    return 0;
  int64_t result = 0;
  for (int64_t k = 0; k * d <= total; ++k) {
    int64_t b1 = math::binomial(n, k);
    int64_t b2 = math::binomial(total - k * d + n - 1, n - 1);
    if (k % 2 == 0)
      result += b1 * b2;
    else
      result -= b1 * b2;
  }
  return result;
}

template <typename bitarray_t>
bitarray_t nth_bounded_partition(int64_t n, int64_t total, int64_t d,
                                 int64_t idx) {
  bitarray_t arr = bits::make_bitarray<bitarray_t>(n);
  int64_t remaining = total;
  for (int64_t slot = n - 1; slot >= 0; --slot) {
    for (int64_t v = 0; v <= std::min(d - 1, remaining); ++v) {
      int64_t cnt = count_bounded_partitions(slot, remaining - v, d);
      if (idx < cnt) {
        arr.set(slot, v);
        remaining -= v;
        break;
      }
      idx -= cnt;
    }
  }
  return arr;
}

template <typename bitarray_t>
int64_t rank_bounded_partition(int64_t n, int64_t total, int64_t d,
                               bitarray_t seq) {
  int64_t result = 0;
  int64_t remaining = total;
  for (int64_t slot = n - 1; slot >= 0; --slot) {
    int64_t v = seq.get(slot);
    for (int64_t vp = 0; vp < v; ++vp)
      result += count_bounded_partitions(slot, remaining - vp, d);
    remaining -= v;
  }
  return result;
}

// clang-format off
#define INSTANTIATE_CBP(BIT_T, NBITS)                                          \
  template bits::BitArray<BIT_T, NBITS>                                        \
  nth_bounded_partition<bits::BitArray<BIT_T, NBITS>>(                         \
      int64_t, int64_t, int64_t, int64_t);                                     \
  template int64_t                                                              \
  rank_bounded_partition<bits::BitArray<BIT_T, NBITS>>(                        \
      int64_t, int64_t, int64_t, bits::BitArray<BIT_T, NBITS>);

#define INSTANTIATE_CBP_ALL(BIT_T)                                             \
  INSTANTIATE_CBP(BIT_T, 1)                                                    \
  INSTANTIATE_CBP(BIT_T, 2)                                                    \
  INSTANTIATE_CBP(BIT_T, 3)                                                    \
  INSTANTIATE_CBP(BIT_T, 4)                                                    \
  INSTANTIATE_CBP(BIT_T, 8)

INSTANTIATE_CBP_ALL(uint16_t)
INSTANTIATE_CBP_ALL(uint32_t)
INSTANTIATE_CBP_ALL(uint64_t)
INSTANTIATE_CBP_ALL(bits::BitsetDynamic)
INSTANTIATE_CBP_ALL(bits::BitsetStatic1)
INSTANTIATE_CBP_ALL(bits::BitsetStatic2)
INSTANTIATE_CBP_ALL(bits::BitsetStatic4)
INSTANTIATE_CBP_ALL(bits::BitsetStatic8)

#undef INSTANTIATE_CBP_ALL
#undef INSTANTIATE_CBP
// clang-format on

} // namespace xdiag::combinatorics
