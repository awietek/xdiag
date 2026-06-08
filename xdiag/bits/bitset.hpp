// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/math/log2.hpp>

namespace xdiag::bits {

// Multi-precision bit storage with arithmetic operations.
//
// Bitset stores a sequence of bits using chunks (uint8_t/16/32/64) and supports
// both dynamic sizing (nchunks=0) and static sizing (nchunks>0). Provides
// bitwise operations, shifts, arithmetic (+,-,*,/), and comparisons, modeling
// unsigned integers.
//
// Template parameters:
//   chunk_tt: Chunk type (uint8_t, uint16_t, uint32_t, uint64_t)
//   nchunks: Number of chunks (0=dynamic using std::vector, >0=static using
//   std::array)
//
// Example:
//   Bitset<uint64_t, 2> bits;      // Static: 128 bits (2 × 64)
//   bits.set(65);                  // Set bit 65 to 1
//   bits = bits + bits;            // Arithmetic operations
//   Bitset<uint64_t, 0> dynamic(200);  // Dynamic: 200 bits
template <typename chunk_tt = uint64_t, int64_t nchunks = 0> class Bitset {
public:
  using chunk_t = chunk_tt;

  // storage_t depends template parameter nchunks:
  // nchunks == 0 -> std::vector<chunk_t>          (dynamic)
  // nchunks != 0 -> std::array<chunk_t, nchunk>   (static)
  using storage_t =
      typename std::conditional<(bool)nchunks, std::array<chunk_t, nchunks>,
                                std::vector<chunk_t>>::type;

  static constexpr size_t nchunkbits = std::numeric_limits<chunk_t>::digits;
  static constexpr size_t nbits = nchunkbits * nchunks;
  static constexpr size_t chunkshift = math::floorlog2(nchunkbits);
  static constexpr chunk_t chunkmask = bitmask<chunk_t>(chunkshift);

  Bitset() = default;
  explicit Bitset(int64_t nbits);

  // Bit-level access — inlined so the compiler sees compile-time constants
  // chunkshift and chunkmask and can strength-reduce the index arithmetic.
  inline bool test(int64_t pos) const noexcept {
    return (chunks_[pos >> chunkshift] >> (pos & chunkmask)) & 1;
  }
  inline void set(int64_t pos) noexcept {
    chunks_[pos >> chunkshift] |= (chunk_t(1) << (pos & chunkmask));
  }
  inline void set(int64_t pos, bool value) noexcept {
    int64_t chunk_idx = pos >> chunkshift;
    int64_t bit_idx = pos & chunkmask;
    if (value) {
      chunks_[chunk_idx] |= (chunk_t(1) << bit_idx);
    } else {
      chunks_[chunk_idx] &= ~(chunk_t(1) << bit_idx);
    }
  }
  inline void reset(int64_t pos) noexcept {
    chunks_[pos >> chunkshift] &= ~(chunk_t(1) << (pos & chunkmask));
  }
  inline void flip(int64_t pos) noexcept {
    chunks_[pos >> chunkshift] ^= (chunk_t(1) << (pos & chunkmask));
  }

  // Bit-level access (ranged) — inlined so the compiler can constant-fold
  // nchunkbits, chunkshift, chunkmask and eliminate branches for typical cases.
  inline void set_range(int64_t start, int64_t length, chunk_t bits) noexcept {
    assert(length <= (int64_t)nchunkbits);
    if (!length) {
      return;
    }
    bits &= bitmask<chunk_t>(length); // mask to length bits
    int64_t end = start + length;
    int64_t startchunk = start >> chunkshift;
    int64_t startbit = start & chunkmask;
    int64_t endchunk = end >> chunkshift;
    int64_t endbit = end & chunkmask;
    if ((endchunk == startchunk) || (endbit == 0)) {
      chunk_t mask = bitmask<chunk_t>(length) << startbit;
      chunks_[startchunk] &= ~mask;
      chunks_[startchunk] |= bits << startbit;
    } else {
      chunk_t negmask1 = bitmask<chunk_t>(startbit);
      chunks_[startchunk] &= negmask1;
      chunks_[startchunk] |= bits << startbit;
      chunk_t mask2 = bitmask<chunk_t>(endbit);
      chunks_[endchunk] &= ~mask2;
      chunks_[endchunk] |= bits >> (nchunkbits - startbit);
    }
  }

  // Atomic OR-into-place variant of set_range. Valid when storage is
  // zero-initialised and each bit position is written by exactly one thread
  // (so the clear step is unnecessary and non-overlapping ORs are safe).
  // Uses reinterpret_cast to std::atomic<chunk_t>*, which is well-defined for
  // lock-free types on all major platforms (GCC, Clang, MSVC / x86-64, ARM64).
  inline void set_range_atomic_or(int64_t start, int64_t length,
                                  chunk_t bits) noexcept {
    static_assert(std::atomic<chunk_t>::is_always_lock_free,
                  "chunk_t must support lock-free atomics for parallel writes");
    bits &= bitmask<chunk_t>(length);
    if (!length)
      return;
    int64_t end = start + length;
    int64_t startchunk = start >> chunkshift;
    int64_t startbit = start & chunkmask;
    int64_t endchunk = end >> chunkshift;
    int64_t endbit = end & chunkmask;
    reinterpret_cast<std::atomic<chunk_t> *>(&chunks_[startchunk])
        ->fetch_or(bits << startbit, std::memory_order_relaxed);
    if (endchunk != startchunk && endbit != 0)
      reinterpret_cast<std::atomic<chunk_t> *>(&chunks_[endchunk])
          ->fetch_or(bits >> (nchunkbits - startbit),
                     std::memory_order_relaxed);
  }

  inline chunk_t get_range(int64_t start, int64_t length) const noexcept {
    assert(length <= (int64_t)nchunkbits);
    if (!length) {
      return chunk_t(0);
    }
    int64_t end = start + length;
    int64_t startchunk = start >> chunkshift;
    int64_t startbit = start & chunkmask;
    int64_t endchunk = end >> chunkshift;
    int64_t endbit = end & chunkmask;
    if ((endchunk == startchunk) || (endbit == 0)) {
      return (chunks_[startchunk] >> startbit) & bitmask<chunk_t>(length);
    } else {
      return ((chunks_[endchunk] & bitmask<chunk_t>(endbit))
              << (nchunkbits - startbit)) |
             (chunks_[startchunk] >> startbit);
    }
  }

  // Bitwise operations
  Bitset operator&(Bitset const &rhs) const;
  Bitset operator|(Bitset const &rhs) const;
  Bitset operator^(Bitset const &rhs) const;
  Bitset operator~() const;

  // Inlined so the compiler sees std::size(chunks_) as a compile-time constant
  // for static arrays (nchunks>0) and can fully unroll the chunk loops.
  inline Bitset &operator&=(Bitset const &rhs) noexcept {
    assert(std::size(chunks_) == std::size(rhs.chunks_));
    for (int64_t i = 0; i < (int64_t)std::size(chunks_); ++i) {
      chunks_[i] &= rhs.chunks_[i];
    }
    return *this;
  }
  inline Bitset &operator|=(Bitset const &rhs) noexcept {
    assert(std::size(chunks_) == std::size(rhs.chunks_));
    for (int64_t i = 0; i < (int64_t)std::size(chunks_); ++i) {
      chunks_[i] |= rhs.chunks_[i];
    }
    return *this;
  }
  inline Bitset &operator^=(Bitset const &rhs) noexcept {
    assert(std::size(chunks_) == std::size(rhs.chunks_));
    for (int64_t i = 0; i < (int64_t)std::size(chunks_); ++i) {
      chunks_[i] ^= rhs.chunks_[i];
    }
    return *this;
  }

  // Shift operations
  Bitset operator<<(int64_t shift) const;
  Bitset operator>>(int64_t shift) const;
  Bitset &operator<<=(int64_t shift) noexcept;
  Bitset &operator>>=(int64_t shift) noexcept;

  // Arithmetic operations
  Bitset operator+(Bitset const &rhs) const;
  Bitset operator-(Bitset const &rhs) const;
  Bitset operator-() const;
  Bitset operator*(Bitset const &rhs) const;
  Bitset operator/(Bitset const &rhs) const;
  Bitset &operator+=(Bitset const &rhs) noexcept;
  Bitset &operator-=(Bitset const &rhs) noexcept;
  Bitset &operator*=(Bitset const &rhs) noexcept;
  Bitset &operator/=(Bitset const &rhs) noexcept;

  // Predicates — inlined for the same reason as operator&=.
  inline bool all() const noexcept {
    for (auto chunk : chunks_) {
      if (chunk != std::numeric_limits<chunk_t>::max()) {
        return false;
      }
    }
    return true;
  }
  inline bool any() const noexcept {
    for (auto chunk : chunks_) {
      if (chunk != 0) {
        return true;
      }
    }
    return false;
  }
  inline bool none() const noexcept { return !any(); }
  int64_t count() const noexcept;

  // Comparisons — inlined so loops over chunks_ unroll for static arrays.
  inline bool operator==(Bitset<chunk_t, nchunks> const &rhs) const noexcept {
    return chunks_ == rhs.chunks_;
  }
  inline bool operator!=(Bitset<chunk_t, nchunks> const &rhs) const noexcept {
    return !operator==(rhs);
  }
  inline bool operator<(Bitset<chunk_t, nchunks> const &rhs) const noexcept {
    assert(std::size(chunks_) == std::size(rhs.chunks_));
    for (int64_t i = (int64_t)std::size(chunks_) - 1; i >= 0; --i) {
      if (chunks_[i] < rhs.chunks_[i]) {
        return true;
      }
      if (chunks_[i] > rhs.chunks_[i]) {
        return false;
      }
    }
    return false; // equal
  }
  inline bool operator<=(Bitset<chunk_t, nchunks> const &rhs) const noexcept {
    return !operator>(rhs);
  }
  inline bool operator>(Bitset<chunk_t, nchunks> const &rhs) const noexcept {
    return rhs.operator<(*this);
  }
  inline bool operator>=(Bitset<chunk_t, nchunks> const &rhs) const noexcept {
    return !operator<(rhs);
  }

  storage_t const &chunks() const noexcept;

private:
  // Optimized shift-by-1 for division algorithm
  void shift_left_by_1() noexcept;
  storage_t chunks_;
};

template <typename chunk_t, int64_t nchunks>
std::string to_string(Bitset<chunk_t, nchunks> const &bits,
                      int64_t size = Bitset<chunk_t, nchunks>::nbits,
                      bool reverse = true);

template <typename chunk_t, int64_t nchunks>
std::ostream &operator<<(std::ostream &out,
                         Bitset<chunk_t, nchunks> const &bits);

// Conversion functions between Bitset and uint64_t (for testing/interop)
template <typename chunk_t = uint64_t, int64_t nchunks = 1>
Bitset<chunk_t, nchunks> make_bitset(uint64_t value);

template <typename chunk_t, int64_t nchunks>
uint64_t to_uint64(Bitset<chunk_t, nchunks> const &bits);

using BitsetDynamic = Bitset<uint64_t, 0>;
using BitsetStatic1 = Bitset<uint64_t, 1>;
using BitsetStatic2 = Bitset<uint64_t, 2>;
using BitsetStatic4 = Bitset<uint64_t, 4>;
using BitsetStatic8 = Bitset<uint64_t, 8>;

} // namespace xdiag::bits

// Specialization of numeric_limits for Bitset (needed for BitArray)
namespace std {
template <> class numeric_limits<xdiag::bits::Bitset<uint64_t, 0>> {
public:
  static constexpr int digits = numeric_limits<int>::max();
};
template <> class numeric_limits<xdiag::bits::Bitset<uint64_t, 1>> {
public:
  static constexpr int digits = numeric_limits<uint64_t>::digits;
};
template <> class numeric_limits<xdiag::bits::Bitset<uint64_t, 2>> {
public:
  static constexpr int digits = numeric_limits<uint64_t>::digits * 2;
};
template <> class numeric_limits<xdiag::bits::Bitset<uint64_t, 4>> {
public:
  static constexpr int digits = numeric_limits<uint64_t>::digits * 4;
};
template <> class numeric_limits<xdiag::bits::Bitset<uint64_t, 8>> {
public:
  static constexpr int digits = numeric_limits<uint64_t>::digits * 8;
};

} // namespace std
