#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

namespace hydra::combinatorics {

template <class int_t> inline int_t get_next_pattern(const int_t &v) {
  // Bit twiddling Hack from
  // http://graphics.stanford.edu/~seander/bithacks.html
  // #NextBitPermutation

  // Fast version (needs __builtin_ctz(v)) (some problem with 0)
  int_t t = v | (v - 1); // t gets v's least significant 0
  return ((t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1)));

  // // Fast version (needs __builtin_ctz(v)) (some problem with 0)
  // int_t t = v | (v - 1); // t gets v's least significant 0
  // return v == 0 ? ~v :((t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) +
  // 1)));

  // // Slow version (should work everywhere)
  // int_t t = (v | (v - 1)) + 1;
  // return v == 0 ? ~v : t | ((((t & -t) / (v & -v)) >> 1) - 1);
}

uint64_t get_nth_pattern(uint64_t n, int n_sites, int n_upspins);
uint64_t get_n_for_pattern(uint64_t pattern, int n_sites, int n_upspins);

} // namespace hydra::combinatorics
