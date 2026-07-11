// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "bitset.hpp"
#include <algorithm>
#include <cassert>
#include <xdiag/bits/popcount.hpp>

// test, set, reset, flip, get_range, set_range,
// operator&=, |=, ^=, all, any, none,
// operator==, !=, <, <=, >, >= are defined inline in bitset.hpp.

namespace xdiag::bits {

template <typename chunk_t>
static constexpr int64_t n_chunks_for_bits(int64_t nbits) {
  constexpr size_t nchunkbits = std::numeric_limits<chunk_t>::digits;
  constexpr size_t chunkshift = math::floorlog2(nchunkbits);
  return nbits > 0 ? ((nbits - 1) >> chunkshift) + 1 : 0;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>::Bitset(int64_t nbits)
    : chunks_([&]() {
        if constexpr (nchunks == 0) {
          // Dynamic: construct vector with calculated size (value-initialized
          // to 0)
          return storage_t(n_chunks_for_bits<chunk_t>(nbits));
        } else {
          // Static: value-initialize array (all elements to 0)
          return storage_t{};
        }
      }()) {
  // For static storage, verify nbits fits in nchunks
  if constexpr (nchunks > 0) {
    assert(n_chunks_for_bits<chunk_t>(nbits) <= nchunks);
  }
}

template <typename chunk_t, int64_t nchunks>
typename Bitset<chunk_t, nchunks>::storage_t const &
Bitset<chunk_t, nchunks>::chunks() const noexcept {
  return chunks_;
}

template <typename chunk_t, int64_t nchunks>
int64_t Bitset<chunk_t, nchunks>::count() const noexcept {
  int64_t total = 0;
  for (auto chunk : chunks_) {
    total += popcount(chunk);
  }
  return total;
}

// Bitwise operations
template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator&(Bitset const &rhs) const {
  Bitset result = *this;
  result &= rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator|(Bitset const &rhs) const {
  Bitset result = *this;
  result |= rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator^(Bitset const &rhs) const {
  Bitset result = *this;
  result ^= rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> Bitset<chunk_t, nchunks>::operator~() const {
  Bitset result = *this;
  for (int64_t i = 0; i < std::size(result.chunks_); ++i) {
    result.chunks_[i] = ~result.chunks_[i];
  }
  return result;
}

// Optimized shift-by-1 for division algorithm (private helper)
template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::shift_left_by_1() noexcept {
  int64_t size = std::size(chunks_);
  chunk_t carry = 0;
  for (int64_t i = 0; i < size; ++i) {
    chunk_t next_carry = chunks_[i] >> (nchunkbits - 1);
    chunks_[i] = (chunks_[i] << 1) | carry;
    carry = next_carry;
  }
}

// Shift operations
template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator<<(int64_t shift) const {
  Bitset result = *this;
  result <<= shift;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator>>(int64_t shift) const {
  Bitset result = *this;
  result >>= shift;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator<<=(int64_t shift) noexcept {
  if (shift == 0)
    return *this;

  int64_t size = std::size(chunks_);
  int64_t chunk_shift = shift >> chunkshift;
  int64_t bit_shift = shift & chunkmask;

  if (chunk_shift >= size) {
    for (auto &chunk : chunks_) {
      chunk = 0;
    }
    return *this;
  }

  if (bit_shift == 0) {
    // Simple chunk-aligned shift
    for (int64_t i = size - 1; i >= chunk_shift; --i) {
      chunks_[i] = chunks_[i - chunk_shift];
    }
    for (int64_t i = 0; i < chunk_shift; ++i) {
      chunks_[i] = 0;
    }
  } else {
    // Shift with bit offset
    int64_t complement_shift = nchunkbits - bit_shift;
    for (int64_t j = size - 1 - chunk_shift; j > 0; --j) {
      chunks_[j + chunk_shift] =
          (chunks_[j] << bit_shift) | (chunks_[j - 1] >> complement_shift);
    }
    chunks_[chunk_shift] = chunks_[0] << bit_shift;
    for (int64_t i = 0; i < chunk_shift; ++i) {
      chunks_[i] = 0;
    }
  }
  return *this;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator>>=(int64_t shift) noexcept {
  if (shift == 0)
    return *this;

  int64_t size = std::size(chunks_);
  int64_t chunk_shift = shift >> chunkshift;
  int64_t bit_shift = shift & chunkmask;

  if (chunk_shift >= size) {
    for (auto &chunk : chunks_) {
      chunk = 0;
    }
    return *this;
  }

  if (bit_shift == 0) {
    // Simple chunk-aligned shift
    for (int64_t i = 0; i < size - chunk_shift; ++i) {
      chunks_[i] = chunks_[i + chunk_shift];
    }
    for (int64_t i = size - chunk_shift; i < size; ++i) {
      chunks_[i] = 0;
    }
  } else {
    // Shift with bit offset
    int64_t complement_shift = nchunkbits - bit_shift;
    for (int64_t i = 0; i < size - chunk_shift - 1; ++i) {
      chunks_[i] = (chunks_[i + chunk_shift] >> bit_shift) |
                   (chunks_[i + chunk_shift + 1] << complement_shift);
    }
    chunks_[size - chunk_shift - 1] = chunks_[size - 1] >> bit_shift;
    for (int64_t i = size - chunk_shift; i < size; ++i) {
      chunks_[i] = 0;
    }
  }
  return *this;
}

// Arithmetic operations
template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator+(Bitset const &rhs) const {
  Bitset result = *this;
  result += rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator-(Bitset const &rhs) const {
  Bitset result = *this;
  result -= rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> Bitset<chunk_t, nchunks>::operator-() const {
  Bitset result = *this;
  // Two's complement: flip all bits and add 1 (correct for unsigned modular
  // arithmetic)
  for (int64_t i = 0; i < std::size(result.chunks_); ++i) {
    result.chunks_[i] = ~result.chunks_[i];
  }
  // Add 1
  chunk_t carry = 1;
  for (int64_t i = 0; i < std::size(result.chunks_) && carry; ++i) {
    chunk_t old_val = result.chunks_[i];
    result.chunks_[i] = old_val + carry;
    carry = (result.chunks_[i] < old_val) ? 1 : 0;
  }
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator*(Bitset const &rhs) const {
  Bitset result = *this;
  result *= rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator/(Bitset const &rhs) const {
  Bitset result = *this;
  result /= rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator+=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));
  chunk_t carry = 0;
  for (int64_t i = 0; i < std::size(chunks_); ++i) {
    chunk_t old_chunk = chunks_[i];
    chunks_[i] += rhs.chunks_[i] + carry;
    // Detect overflow: if result is less than either operand, we had overflow
    carry =
        (chunks_[i] < old_chunk || (carry && chunks_[i] == old_chunk)) ? 1 : 0;
  }
  return *this;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator-=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));
  chunk_t borrow = 0;
  for (int64_t i = 0; i < std::size(chunks_); ++i) {
    chunk_t old_chunk = chunks_[i];
    chunks_[i] -= rhs.chunks_[i] + borrow;
    // Detect underflow: if result is greater than original, we had underflow
    borrow =
        (chunks_[i] > old_chunk || (borrow && chunks_[i] == old_chunk)) ? 1 : 0;
  }
  return *this;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator*=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));
  Bitset<chunk_t, nchunks> result(std::size(chunks_) * nchunkbits);

  // Grade-school multiplication without extended types
  // Split each chunk multiplication into high and low parts
  constexpr int half_bits = nchunkbits / 2;
  constexpr chunk_t low_mask = (static_cast<chunk_t>(1) << half_bits) - 1;

  for (int64_t i = 0; i < std::size(chunks_); ++i) {
    if (rhs.chunks_[i] == 0)
      continue;

    chunk_t a_lo = rhs.chunks_[i] & low_mask;
    chunk_t a_hi = rhs.chunks_[i] >> half_bits;

    chunk_t carry = 0;
    for (int64_t j = 0; j < std::size(chunks_) - i; ++j) {
      chunk_t b_lo = chunks_[j] & low_mask;
      chunk_t b_hi = chunks_[j] >> half_bits;

      // Compute a * b as (a_hi * 2^k + a_lo) * (b_hi * 2^k + b_lo)
      chunk_t p_ll = a_lo * b_lo;
      chunk_t p_lh = a_lo * b_hi;
      chunk_t p_hl = a_hi * b_lo;
      chunk_t p_hh = a_hi * b_hi;

      // Combine: p_ll + (p_lh + p_hl) * 2^k + p_hh * 2^(2k)
      chunk_t mid = p_lh + p_hl;
      chunk_t mid_carry = (mid < p_lh) ? 1 : 0;

      chunk_t low = p_ll + ((mid & low_mask) << half_bits);
      chunk_t low_carry = (low < p_ll) ? 1 : 0;

      chunk_t high =
          p_hh + (mid >> half_bits) + (mid_carry << half_bits) + low_carry;

      // Add to result with carry
      chunk_t old_val = result.chunks_[i + j];
      chunk_t temp = old_val + low;
      chunk_t add_carry1 = (temp < old_val) ? 1 : 0;
      result.chunks_[i + j] = temp + carry;
      chunk_t add_carry2 = (result.chunks_[i + j] < temp) ? 1 : 0;
      carry = high + add_carry1 + add_carry2;
    }
  }

  *this = result;
  return *this;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator/=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));

  // Division by zero check
  if (rhs.none()) {
    // Division by zero - return 0
    for (auto &chunk : chunks_) {
      chunk = 0;
    }
    return *this;
  }

  // Long division algorithm for multi-precision integers
  Bitset<chunk_t, nchunks> quotient(std::size(chunks_) * nchunkbits);
  Bitset<chunk_t, nchunks> remainder(std::size(chunks_) * nchunkbits);

  // Find the most significant bit position
  int64_t nbits = std::size(chunks_) * nchunkbits;

  for (int64_t i = nbits - 1; i >= 0; --i) {
    // Shift remainder left by 1 (optimized)
    remainder.shift_left_by_1();
    // Set the least significant bit of remainder to bit i of dividend
    if (test(i)) {
      remainder.set(0);
    }

    // If remainder >= divisor, subtract divisor and set quotient bit
    if (remainder >= rhs) {
      remainder -= rhs;
      quotient.set(i);
    }
  }

  *this = quotient;
  return *this;
}

template <typename chunk_t, int64_t nchunks>
std::string to_string(Bitset<chunk_t, nchunks> const &bits, int64_t size,
                      bool reverse) {
  std::string str;
  for (int64_t i = 0; i < size; ++i) {
    str += bits.test(i) ? std::string("1") : std::string("0");
  }
  if (reverse) {
    std::reverse(str.begin(), str.end());
  }
  return str;
}

template <typename chunk_t, int64_t nchunks>
std::ostream &operator<<(std::ostream &out,
                         Bitset<chunk_t, nchunks> const &bits) {
  out << to_string(bits);
  return out;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> make_bitset(uint64_t value) {
  Bitset<chunk_t, nchunks> bits(64);
  for (int i = 0; i < 64; ++i) {
    if (value & (1ULL << i)) {
      bits.set(i);
    }
  }
  return bits;
}

template <typename chunk_t, int64_t nchunks>
uint64_t to_uint64(Bitset<chunk_t, nchunks> const &bits) {
  uint64_t result = 0;
  constexpr int64_t chunk_bits = std::numeric_limits<chunk_t>::digits;

  int64_t max_bits;
  if constexpr (nchunks == 0) {
    // Dynamic: runtime calculation
    max_bits =
        std::min(int64_t(64), int64_t(std::size(bits.chunks()) * chunk_bits));
  } else {
    // Static: compile-time calculation
    constexpr int64_t max_bits_static =
        std::min(int64_t(64), nchunks * chunk_bits);
    max_bits = max_bits_static;
  }

  for (int i = 0; i < max_bits; ++i) {
    if (bits.test(i)) {
      result |= (1ULL << i);
    }
  }
  return result;
}

} // namespace xdiag::bits

using namespace xdiag::bits;

// Explicit template instantiations
#define INSTANTIATE_XDIAG_BITS_BITSET(CHUNK_T, NCHUNKS)                        \
  template class xdiag::bits::Bitset<CHUNK_T, NCHUNKS>;                        \
  template std::string xdiag::bits::to_string(                                 \
      Bitset<CHUNK_T, NCHUNKS> const &, int64_t, bool);                        \
  template std::ostream &xdiag::bits::operator<<(                              \
      std::ostream &, Bitset<CHUNK_T, NCHUNKS> const &);                       \
  template Bitset<CHUNK_T, NCHUNKS> xdiag::bits::make_bitset(uint64_t);        \
  template uint64_t xdiag::bits::to_uint64(Bitset<CHUNK_T, NCHUNKS> const &);


// BEGIN_INSTANTIATION_GROUP(uint32_t)
INSTANTIATE_XDIAG_BITS_BITSET(uint32_t, 0)
INSTANTIATE_XDIAG_BITS_BITSET(uint32_t, 1)
INSTANTIATE_XDIAG_BITS_BITSET(uint32_t, 2)
INSTANTIATE_XDIAG_BITS_BITSET(uint32_t, 4)
INSTANTIATE_XDIAG_BITS_BITSET(uint32_t, 8)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(uint64_t)
INSTANTIATE_XDIAG_BITS_BITSET(uint64_t, 0)
INSTANTIATE_XDIAG_BITS_BITSET(uint64_t, 1)
INSTANTIATE_XDIAG_BITS_BITSET(uint64_t, 2)
INSTANTIATE_XDIAG_BITS_BITSET(uint64_t, 4)
INSTANTIATE_XDIAG_BITS_BITSET(uint64_t, 8)
// END_INSTANTIATION_GROUP

#undef INSTANTIATE_XDIAG_BITS_BITSET
