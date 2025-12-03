// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "bitset.hpp"
#include <bitset>
#include <cassert>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::bits {

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>::Bitset(uint64_t bits) {
  int64_t nbits = std::numeric_limits<uint64_t>::digits;
  int64_t size = (std::max(nbits - 1, (int64_t)0)) / nchunkbits + 1;
  resize(size);
  uint64_t mask = ~(chunk_t)0;
  int64_t nchunk = 0;
  while ((bool)bits && (nchunk < std::size(chunks_))) {
    chunks_[nchunk++] = bits & mask;
    if constexpr (nchunkbits == std::numeric_limits<uint64_t>::digits) {
      break;
    } else {
      bits >>= nchunkbits;
    }
  }
}

template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::resize(int64_t n_chunks) {
  if constexpr (nchunks == 0) {
    chunks_.resize(n_chunks);
  }
}

template <typename chunk_t, int64_t nchunks>
int64_t Bitset<chunk_t, nchunks>::size() const {
  return std::size(chunks_);
}
template <typename chunk_t, int64_t nchunks>
typename Bitset<chunk_t, nchunks>::iterator_t
Bitset<chunk_t, nchunks>::begin() const {
  return std::begin(chunks_);
}
template <typename chunk_t, int64_t nchunks>
typename Bitset<chunk_t, nchunks>::iterator_t
Bitset<chunk_t, nchunks>::end() const {
  return std::end(chunks_);
}

template <typename chunk_t, int64_t nchunks>
chunk_t Bitset<chunk_t, nchunks>::get_bits(int64_t start,
                                           int64_t length) const {
  assert(length <= nchunkbits);
  int64_t end = start + length;
  int64_t startchunk = start >> chunkshift; // divide by nchunkbits
  int64_t startbit = start & chunkmask;     // modulo    nchunkbits
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

template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::set_bits(int64_t start, int64_t length,
                                        chunk_t bits) {
  assert(length <= nchunkbits);
  int64_t end = start + length;
  int64_t startchunk = start >> chunkshift; // divide by nchunkbits
  int64_t startbit = start & chunkmask;     // modulo    nchunkbits
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

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::operator==(
    Bitset<chunk_t, nchunks> const &rhs) const {
  return chunks_ == rhs.chunks_;
}

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::operator!=(
    Bitset<chunk_t, nchunks> const &rhs) const {
  return !operator==(rhs);
}

std::string to_string(uint64_t bits) {
  constexpr size_t nbits = std::numeric_limits<uint64_t>::digits;
  return std::bitset<nbits>(bits).to_string();
}

template <typename chunk_t, int64_t nchunks>
std::string to_string(Bitset<chunk_t, nchunks> const &bits) {
  std::string str;
  for (auto chunk : bits) {
    str = std::bitset<Bitset<chunk_t, nchunks>::nchunkbits>(chunk).to_string() +
          str;
  }
  return str;
}

template <typename chunk_t, int64_t nchunks>
std::ostream &operator<<(std::ostream &out,
                         Bitset<chunk_t, nchunks> const &bits) {
  out << to_string(bits);
  return out;
}

template class Bitset<uint64_t, 0>;
template class Bitset<uint64_t, 1>;
template class Bitset<uint64_t, 2>;
template class Bitset<uint64_t, 3>;
template class Bitset<uint64_t, 4>;

template class Bitset<uint32_t, 0>;
template class Bitset<uint32_t, 1>;
template class Bitset<uint32_t, 2>;
template class Bitset<uint32_t, 3>;
template class Bitset<uint32_t, 4>;

template class Bitset<uint16_t, 0>;
template class Bitset<uint16_t, 1>;
template class Bitset<uint16_t, 2>;
template class Bitset<uint16_t, 3>;
template class Bitset<uint16_t, 4>;

template class Bitset<uint8_t, 0>;
template class Bitset<uint8_t, 1>;
template class Bitset<uint8_t, 2>;
template class Bitset<uint8_t, 3>;
template class Bitset<uint8_t, 4>;

template std::string to_string(Bitset<uint64_t, 0> const &);
template std::string to_string(Bitset<uint64_t, 1> const &);
template std::string to_string(Bitset<uint64_t, 2> const &);
template std::string to_string(Bitset<uint64_t, 3> const &);
template std::string to_string(Bitset<uint64_t, 4> const &);

template std::string to_string(Bitset<uint32_t, 0> const &);
template std::string to_string(Bitset<uint32_t, 1> const &);
template std::string to_string(Bitset<uint32_t, 2> const &);
template std::string to_string(Bitset<uint32_t, 3> const &);
template std::string to_string(Bitset<uint32_t, 4> const &);

template std::string to_string(Bitset<uint16_t, 0> const &);
template std::string to_string(Bitset<uint16_t, 1> const &);
template std::string to_string(Bitset<uint16_t, 2> const &);
template std::string to_string(Bitset<uint16_t, 3> const &);
template std::string to_string(Bitset<uint16_t, 4> const &);

template std::string to_string(Bitset<uint8_t, 0> const &);
template std::string to_string(Bitset<uint8_t, 1> const &);
template std::string to_string(Bitset<uint8_t, 2> const &);
template std::string to_string(Bitset<uint8_t, 3> const &);
template std::string to_string(Bitset<uint8_t, 4> const &);

template std::ostream &operator<<(std::ostream &, Bitset<uint64_t, 0> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint64_t, 1> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint64_t, 2> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint64_t, 3> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint64_t, 4> const &);

template std::ostream &operator<<(std::ostream &, Bitset<uint32_t, 0> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint32_t, 1> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint32_t, 2> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint32_t, 3> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint32_t, 4> const &);

template std::ostream &operator<<(std::ostream &, Bitset<uint16_t, 0> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint16_t, 1> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint16_t, 2> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint16_t, 3> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint16_t, 4> const &);

template std::ostream &operator<<(std::ostream &, Bitset<uint8_t, 0> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint8_t, 1> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint8_t, 2> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint8_t, 3> const &);
template std::ostream &operator<<(std::ostream &, Bitset<uint8_t, 4> const &);
} // namespace xdiag::bits
