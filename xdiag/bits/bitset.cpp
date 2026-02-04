// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "bitset.hpp"
#include <bitset>
#include <cassert>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/gbit.hpp>
#include <xdiag/bits/popcnt.hpp>
#include <xdiag/utils/logger.hpp>

#include <xdiag/bits/log2.hpp>

namespace xdiag::bits {

template <typename chunk_t>
static constexpr int64_t n_chunks_for_bits(int64_t nbits) {
  constexpr size_t nchunkbits = std::numeric_limits<chunk_t>::digits;
  constexpr size_t chunkshift = floorlog2(nchunkbits);
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

// Bit-level access
template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::test(int64_t pos) const noexcept {
  int64_t chunk_idx = pos >> chunkshift_;
  int64_t bit_idx = pos & chunkmask_;
  return gbit(chunks_[chunk_idx], bit_idx);
}

template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::set(int64_t pos) noexcept {
  int64_t chunk_idx = pos >> chunkshift_;
  int64_t bit_idx = pos & chunkmask_;
  chunks_[chunk_idx] |= ((chunk_t)1 << bit_idx);
}

template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::set(int64_t pos, bool value) noexcept {
  int64_t chunk_idx = pos >> chunkshift_;
  int64_t bit_idx = pos & chunkmask_;
  if (value) {
    chunks_[chunk_idx] |= ((chunk_t)1 << bit_idx);
  } else {
    chunks_[chunk_idx] &= ~((chunk_t)1 << bit_idx);
  }
}

template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::reset(int64_t pos) noexcept {
  int64_t chunk_idx = pos >> chunkshift_;
  int64_t bit_idx = pos & chunkmask_;
  chunks_[chunk_idx] &= ~((chunk_t)1 << bit_idx);
}

template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::flip(int64_t pos) noexcept {
  int64_t chunk_idx = pos >> chunkshift_;
  int64_t bit_idx = pos & chunkmask_;
  chunks_[chunk_idx] ^= ((chunk_t)1 << bit_idx);
}

// Bit-level access (ranged)
template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::set_range(int64_t start, int64_t length,
                                         chunk_t bits) noexcept {
  assert(length <= nchunkbits_);
  if (!length) {
    return;
  }
  int64_t end = start + length;
  int64_t startchunk = start >> chunkshift_; // divide by nchunkbits
  int64_t startbit = start & chunkmask_;     // modulo    nchunkbits
  int64_t endchunk = end >> chunkshift_;
  int64_t endbit = end & chunkmask_;
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
    chunks_[endchunk] |= bits >> (nchunkbits_ - startbit);
  }
}

template <typename chunk_t, int64_t nchunks>
chunk_t Bitset<chunk_t, nchunks>::get_range(int64_t start,
                                            int64_t length) const noexcept {
  assert(length <= nchunkbits_);
  if (!length) {
    return (chunk_t)0;
  }
  int64_t end = start + length;
  int64_t startchunk = start >> chunkshift_; // divide by nchunkbits
  int64_t startbit = start & chunkmask_;     // modulo    nchunkbits
  int64_t endchunk = end >> chunkshift_;
  int64_t endbit = end & chunkmask_;
  if ((endchunk == startchunk) || (endbit == 0)) {
    return (chunks_[startchunk] >> startbit) & bitmask<chunk_t>(length);
  } else {
    return ((chunks_[endchunk] & bitmask<chunk_t>(endbit))
            << (nchunkbits_ - startbit)) |
           (chunks_[startchunk] >> startbit);
  }
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

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator&=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));
  for (int64_t i = 0; i < std::size(chunks_); ++i) {
    chunks_[i] &= rhs.chunks_[i];
  }
  return *this;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator|=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));
  for (int64_t i = 0; i < std::size(chunks_); ++i) {
    chunks_[i] |= rhs.chunks_[i];
  }
  return *this;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator^=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));
  for (int64_t i = 0; i < std::size(chunks_); ++i) {
    chunks_[i] ^= rhs.chunks_[i];
  }
  return *this;
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
  int64_t chunk_shift = shift >> chunkshift_;
  int64_t bit_shift = shift & chunkmask_;

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
    int64_t complement_shift = nchunkbits_ - bit_shift;
    for (int64_t i = size - 1; i > chunk_shift; --i) {
      chunks_[i] = (chunks_[i - chunk_shift] << bit_shift) |
                   (chunks_[i - chunk_shift - 1] >> complement_shift);
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
  int64_t chunk_shift = shift >> chunkshift_;
  int64_t bit_shift = shift & chunkmask_;

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
    int64_t complement_shift = nchunkbits_ - bit_shift;
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

// Predicates
template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::all() const noexcept {
  for (auto chunk : chunks_) {
    if (chunk != std::numeric_limits<chunk_t>::max()) {
      return false;
    }
  }
  return true;
}

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::any() const noexcept {
  for (auto chunk : chunks_) {
    if (chunk != 0) {
      return true;
    }
  }
  return false;
}

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::none() const noexcept {
  return !any();
}

template <typename chunk_t, int64_t nchunks>
int64_t Bitset<chunk_t, nchunks>::count() const noexcept {
  int64_t total = 0;
  for (auto chunk : chunks_) {
    total += popcnt(chunk);
  }
  return total;
}

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::operator==(
    Bitset<chunk_t, nchunks> const &rhs) const noexcept {
  return chunks_ == rhs.chunks_;
}

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::operator!=(
    Bitset<chunk_t, nchunks> const &rhs) const noexcept {
  return !operator==(rhs);
}

template <typename chunk_t, int64_t nchunks>
std::string to_string(Bitset<chunk_t, nchunks> const &bits) {
  std::string str;
  for (auto const &chunk : bits.chunks()) {
    str = std::bitset<bits.nchunkbits()>(chunk).to_string() + str;
  }
  return str;
}

template <typename chunk_t, int64_t nchunks>
std::ostream &operator<<(std::ostream &out,
                         Bitset<chunk_t, nchunks> const &bits) {
  out << to_string(bits);
  return out;
}

// Explicit template instantiations
#define INSTANTIATE_XDIAG_BITS_BITSET(CHUNK_T, NCHUNKS)                        \
  template class Bitset<CHUNK_T, NCHUNKS>;                                     \
  template std::string to_string(Bitset<CHUNK_T, NCHUNKS> const &);            \
  template std::ostream &operator<<(std::ostream &,                            \
                                    Bitset<CHUNK_T, NCHUNKS> const &);

#define INSTANTIATE_XDIAG_BITS_BITSET_FOR_NCHUNKS(CHUNK_T)                     \
  INSTANTIATE_XDIAG_BITS_BITSET(CHUNK_T, 0)                                    \
  INSTANTIATE_XDIAG_BITS_BITSET(CHUNK_T, 1)                                    \
  INSTANTIATE_XDIAG_BITS_BITSET(CHUNK_T, 2)                                    \
  INSTANTIATE_XDIAG_BITS_BITSET(CHUNK_T, 4)                                    \
  INSTANTIATE_XDIAG_BITS_BITSET(CHUNK_T, 8)

INSTANTIATE_XDIAG_BITS_BITSET_FOR_NCHUNKS(uint8_t)
INSTANTIATE_XDIAG_BITS_BITSET_FOR_NCHUNKS(uint16_t)
INSTANTIATE_XDIAG_BITS_BITSET_FOR_NCHUNKS(uint32_t)
INSTANTIATE_XDIAG_BITS_BITSET_FOR_NCHUNKS(uint64_t)

#undef INSTANTIATE_XDIAG_BITS_BITSET_FOR_NCHUNKS
#undef INSTANTIATE_XDIAG_BITS_BITSET

} // namespace xdiag::bits
