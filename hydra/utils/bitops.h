#ifndef HYDRA_UTILS_BITOPS_
#define HYDRA_UTILS_BITOPS_

#include <hydra/common.h>

namespace hydra {
namespace utils {

// gbit
template <class uint_t> inline uint_t gbit(uint_t x, int n) {
  return (x >> n) & (uint_t)1;
}
template <> inline uint16 gbit<uint16>(uint16 x, int n) {
  return (x >> n) & (1U);
}
template <> inline uint32 gbit<uint32>(uint32 x, int n) {
  return (x >> n) & (1U);
}
template <> inline uint64 gbit<uint64>(uint64 x, int n) {
  return (x >> n) & (1UL);
}

// gbits
template <class uint_t> inline uint_t gbits(uint_t x, int m, int n) {
  return (x >> n) & (((uint_t)1 << m) - 1);
}
template <> inline uint16 gbits<uint16>(uint16 x, int m, int n) {
  return (x >> n) & ((1U << m) - 1);
}
template <> inline uint32 gbits<uint32>(uint32 x, int m, int n) {
  return (x >> n) & ((1U << m) - 1);
}
template <> inline uint64 gbits<uint64>(uint64 x, int m, int n) {
  return (x >> n) & ((1UL << m) - 1);
}

// sbit
template <class uint_t> inline void sbit(uint_t &x, int n, uint_t b) {
  x &= ~(1 << n) | (b << n);
}

// popcnt
template <class uint_t> inline int popcnt(uint_t x) {
  return __builtin_popcountll(x);
}
template <> inline int popcnt<uint16>(uint16 x) {
  return __builtin_popcount(x);
}
template <> inline int popcnt<uint32>(uint32 x) {
  return __builtin_popcount(x);
}
template <> inline int popcnt<uint64>(uint64 x) {
  return __builtin_popcountll(x);
}

} // namespace utils
} // namespace hydra
#endif
