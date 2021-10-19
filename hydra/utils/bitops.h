#pragma once

#include <hydra/common.h>
#include <sstream>
#include <string>

namespace hydra::utils {

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

// bits_to_string
template <typename bit_t>
std::string bits_to_string(bit_t bits, int n, bool reverse = true) {
  std::stringstream s;
  for (int i = 0; i < n; ++i)
    s << utils::gbit(bits, i);
  std::string st = s.str();
  return reverse ? std::string(st.rbegin(), st.rend()) : st;
}

} // namespace hydra::utils
