// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>

namespace xdiag::bits {

// get
template <class bit_t> constexpr int64_t get(bit_t x, uint32_t n) {
  return (x >> n) & (bit_t)1;
}

template <class chunk_t, int64_t nchunks>
inline int64_t get(Bitset<chunk_t, nchunks> const &bits, int64_t n) {
  return bits.test(n) ? 1 : 0;
}

template <class storage_t, int nbits>
inline int64_t get(BitArray<storage_t, nbits> const &arr, int64_t n) {
  return arr.get(n);
}

// get_range
template <class bit_t>
constexpr bit_t get_range(bit_t src, uint32_t start, uint32_t len) {
  return (src >> start) & (((bit_t)1 << len) - 1);
}

template <class chunk_t, int64_t nchunks>
inline chunk_t get_range(Bitset<chunk_t, nchunks> const &bits, int64_t start,
                         int64_t len) {
  return bits.get_range(start, len);
}

// set (2-arg: set bit/slot to 1)
template <class bit_t> constexpr void set(bit_t &x, int64_t b) {
  x |= bit_t(1) << b;
}
template <class chunk_t, int64_t nchunks>
inline void set(Bitset<chunk_t, nchunks> &bits, int64_t n) {
  bits.set(n);
}
template <class storage_t, int nbits>
inline void set(BitArray<storage_t, nbits> &arr, int64_t n) {
  arr.set(n, 1);
}

// set (3-arg: set bit/slot to val)
template <class bit_t> constexpr void set(bit_t &x, int64_t b, int64_t val) {
  bit_t mask = bit_t(1) << b;
  x = (x & ~mask) | (bit_t(-(val & 1)) & mask);
}
template <class chunk_t, int64_t nchunks>
inline void set(Bitset<chunk_t, nchunks> &bits, int64_t n, int64_t val) {
  bits.set(n, (bool)(val & 1));
}
template <class storage_t, int nbits>
inline void set(BitArray<storage_t, nbits> &arr, int64_t n, int64_t val) {
  arr.set(n, val);
}

} // namespace xdiag::bits
