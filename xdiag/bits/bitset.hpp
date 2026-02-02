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

#include <xdiag/bits/log2.hpp>
#include <xdiag/common.hpp>

namespace xdiag::bits {

template <typename chunk_t = uint64_t, int64_t nchunks = 0> class Bitset {
public:
  static constexpr size_t nchunkbits = std::numeric_limits<chunk_t>::digits;
  static constexpr size_t chunkshift = floorlog2(nchunkbits);
  static constexpr chunk_t chunkmask = ((chunk_t)1 << chunkshift) - 1;
  using storage_t =
      typename std::conditional<(bool)nchunks, std::array<chunk_t, nchunks>,
                                std::vector<chunk_t>>::type;
  using iterator_t = typename storage_t::const_iterator;

  Bitset() = default;
  Bitset(uint64_t bits);
  void resize(int64_t n_chunks); // used when nchunks=0 -> vector

  int64_t size() const;     // returns number of chunks
  iterator_t begin() const; // iteration over chunks
  iterator_t end() const;

  void set_bits(int64_t start, int64_t length, chunk_t bits);
  chunk_t get_bits(int64_t start, int64_t length) const;

  bool operator==(Bitset<chunk_t, nchunks> const &rhs) const;
  bool operator!=(Bitset<chunk_t, nchunks> const &rhs) const;

private:
  storage_t chunks_;
};

std::string to_string(uint64_t bits);

template <typename chunk_t, int64_t nchunks>
std::string to_string(Bitset<chunk_t, nchunks> const &bits);

template <typename chunk_t, int64_t nchunks>
std::ostream &operator<<(std::ostream &out,
                         Bitset<chunk_t, nchunks> const &bits);

} // namespace xdiag::bits
