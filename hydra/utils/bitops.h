#pragma once

#include <hydra/common.h>
#include <immintrin.h>
#include <sstream>
#include <string>

namespace hydra::bitops {

// gbit
template <class uint_t> inline uint_t gbit(uint_t x, int n) {
  return (x >> n) & (uint_t)1;
}

// gbits
template <class uint_t> inline uint_t gbits(uint_t x, int m, int n) {
  return (x >> n) & (((uint_t)1 << m) - 1);
}

// sbit
template <class uint_t> inline void sbit(uint_t &x, int n, uint_t b) {
  x &= ~((uint_t)1 << n) | (b << n);
}

// popcnt
inline int popcnt(int x) { return __builtin_popcount(x); }
inline int popcnt(uint16 x) { return __builtin_popcount(x); }
inline int popcnt(uint32 x) { return __builtin_popcount(x); }
inline int popcnt(uint64 x) { return __builtin_popcountll(x); }

// extract -> PEXT instruction
inline uint16 extract(uint16 src, uint16 mask) { return _pext_u32(src, mask); }
inline uint32 extract(uint32 src, uint32 mask) { return _pext_u32(src, mask); }
inline uint64 extract(uint64 src, uint64 mask) { return _pext_u64(src, mask); }

// deposit -> PDEP instruction
inline uint16 deposit(uint16 src, uint16 mask) { return _pdep_u32(src, mask); }
inline uint32 deposit(uint32 src, uint32 mask) { return _pdep_u32(src, mask); }
inline uint64 deposit(uint64 src, uint64 mask) { return _pdep_u64(src, mask); }

// bits_to_string
template <typename bit_t>
std::string bits_to_string(bit_t bits, int n, bool reverse = true) {
  std::stringstream s;
  for (int i = 0; i < n; ++i)
    s << gbit(bits, i);
  std::string st = s.str();
  return reverse ? std::string(st.rbegin(), st.rend()) : st;
}

} // namespace hydra::bitops
