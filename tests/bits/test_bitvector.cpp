// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/bitvector.hpp>
#include <xdiag/utils/logger.hpp>

#include <limits>
#include <random>

using namespace xdiag;
using namespace xdiag::bits;

template <typename bit_t> void test_bitvector() {
  std::mt19937 rng(42);

  for (int nbits = 0; nbits <= std::numeric_limits<bit_t>::digits; ++nbits) {
    for (int size = 0; size < 100; ++size) {
      auto vec = BitVector<bit_t>(nbits, size);
      std::uniform_int_distribution<bit_t> dist(0, bitmask<bit_t>(nbits));
      for (int idx = 0; idx < size; ++idx) {
        for (int r = 0; r < 4; ++r) {
          bit_t bits = dist(rng);
          vec[idx] = bits;
          bit_t bits2 = vec[idx];
          REQUIRE(bits == bits2);
        }
      }
    }
  }
}

TEST_CASE("bitvector", "[bits]") try {
  Log("Testing bitvector uint8_t");
  test_bitvector<uint8_t>();
  Log("Testing bitvector uint16_t");
  test_bitvector<uint16_t>();
  Log("Testing bitvector uint32_t");
  test_bitvector<uint32_t>();
  Log("Testing bitvector uint64_t");
  test_bitvector<uint64_t>();
} catch (xdiag::Error const &e) {
  error_trace(e);
}
