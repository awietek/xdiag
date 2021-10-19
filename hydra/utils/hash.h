#pragma once
#include <hydra/common.h>

namespace hydra::utils {

inline uint16 fast_hash(uint16 x) {
  x = ((x >> 8) ^ x) * 0x9f3b;
  x = ((x >> 8) ^ x) * 0x9f3b;
  x = (x >> 8) ^ x;
  return x;
}
inline uint32 fast_hash(uint32 x) {
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = (x >> 16) ^ x;
  return x;
}
inline uint64 fast_hash(uint64 x) {
  x = (x ^ (x >> 30)) * uint64(0xbf58476d1ce4e5b9);
  x = (x ^ (x >> 27)) * uint64(0x94d049bb133111eb);
  x = x ^ (x >> 31);
  return x;
}

} // namespace hydra::utils
