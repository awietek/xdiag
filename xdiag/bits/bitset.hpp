// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/log2.hpp>
#include <xdiag/common.hpp>

namespace xdiag::bits {

template <typename chunk_t = uint64_t, int64_t nchunks = 0> class Bitset {
public:
  // storage_t depends template parameter nchunks:
  // nchunks == 0 -> std::vector<chunk_t>          (dynamic)
  // nchunks != 0 -> std::array<chunk_t, nchunk>   (static)
  using storage_t =
      typename std::conditional<(bool)nchunks, std::array<chunk_t, nchunks>,
                                std::vector<chunk_t>>::type;

  Bitset() = default;
  explicit Bitset(int64_t nbits);

  // Bit-level access
  bool test(int64_t pos) const noexcept;
  void set(int64_t pos) noexcept;             // fast: always sets to 1
  void set(int64_t pos, bool value) noexcept; // conditional set
  void reset(int64_t pos) noexcept;
  void flip(int64_t pos) noexcept;

  // Bit-level access (ranged), assuming length <= nchunkbits_
  void set_range(int64_t start, int64_t length, chunk_t bits) noexcept;
  chunk_t get_range(int64_t start, int64_t length) const noexcept;

  // Bitwise operations
  Bitset operator&(Bitset const &rhs) const;
  Bitset operator|(Bitset const &rhs) const;
  Bitset operator^(Bitset const &rhs) const;
  Bitset operator~() const;
  Bitset &operator&=(Bitset const &rhs) noexcept;
  Bitset &operator|=(Bitset const &rhs) noexcept;
  Bitset &operator^=(Bitset const &rhs) noexcept;

  // Shift operations
  Bitset operator<<(int64_t shift) const;
  Bitset operator>>(int64_t shift) const;
  Bitset &operator<<=(int64_t shift) noexcept;
  Bitset &operator>>=(int64_t shift) noexcept;

  // Predicates
  bool all() const noexcept;
  bool any() const noexcept;
  bool none() const noexcept;
  int64_t count() const noexcept;

  bool operator==(Bitset<chunk_t, nchunks> const &rhs) const noexcept;
  bool operator!=(Bitset<chunk_t, nchunks> const &rhs) const noexcept;

  storage_t const &chunks() const noexcept;
  static constexpr size_t nchunkbits() noexcept { return nchunkbits_; }

private:
  static constexpr size_t nchunkbits_ = std::numeric_limits<chunk_t>::digits;
  static constexpr size_t chunkshift_ = floorlog2(nchunkbits_);
  static constexpr chunk_t chunkmask_ = bitmask<chunk_t>(chunkshift_);
  storage_t chunks_;
};

template <typename chunk_t, int64_t nchunks>
std::string to_string(Bitset<chunk_t, nchunks> const &bits);

template <typename chunk_t, int64_t nchunks>
std::ostream &operator<<(std::ostream &out,
                         Bitset<chunk_t, nchunks> const &bits);

} // namespace xdiag::bits
