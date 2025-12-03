// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/utils/logger.hpp>

#include <limits>
#include <random>

using namespace xdiag;
using namespace xdiag::bits;

template <typename chunk_t> void test_bitset_get_bits() {
  constexpr int chunkdigits = std::numeric_limits<chunk_t>::digits;
  constexpr int nchunks = 64 / chunkdigits;

  // Random uint64_t
  std::mt19937 rng(42);
  std::uniform_int_distribution<uint64_t> dist(
      std::numeric_limits<std::uint64_t>::min(),
      std::numeric_limits<std::uint64_t>::max());

  std::uniform_int_distribution<std::mt19937::result_type> diststart(0, 64);
  std::uniform_int_distribution<std::mt19937::result_type> distlength(
      0, chunkdigits);

  for (int i = 0; i < 500; ++i) {
    uint64_t random_uint64 = dist(rng);
    int64_t start = diststart(rng);
    int64_t length = std::min((int64_t)distlength(rng), 64 - start);
    chunk_t random_uint64_filtered =
        (random_uint64 >> start) & bitmask<chunk_t>(length);

    if constexpr (chunkdigits > 8) {
      auto bits = Bitset<chunk_t, nchunks>(random_uint64);
      chunk_t bits_filtered = bits.get_bits(start, length);
      REQUIRE(random_uint64_filtered == bits_filtered);
    }
    auto bits = Bitset<chunk_t>(random_uint64);
    chunk_t bits_filtered = bits.get_bits(start, length);
    REQUIRE(random_uint64_filtered == bits_filtered);
  }
}

template <typename chunk_t> void test_bitset_set_bits() {
  constexpr int chunkdigits = std::numeric_limits<chunk_t>::digits;
  constexpr int nchunks = 64 / chunkdigits;
  auto bitsstatic = Bitset<chunk_t, nchunks>();
  auto bitsdynamic = Bitset<chunk_t>();
  bitsdynamic.resize(nchunks);

  // Random uint64_t
  std::mt19937 rng(42);

  std::uniform_int_distribution<std::mt19937::result_type> diststart(0, 64);
  std::uniform_int_distribution<std::mt19937::result_type> distlength(
      0, chunkdigits);

  for (int i = 0; i < 100; ++i) {
    int64_t start = diststart(rng);
    int64_t length = std::min((int64_t)distlength(rng), 64 - start);
    std::uniform_int_distribution<chunk_t> dist(0, ((chunk_t)1 << length) - 1);
    chunk_t random_bits = dist(rng);
    if constexpr (chunkdigits > 8) {
      bitsstatic.set_bits(start, length, random_bits);
      chunk_t random_bits2 = bitsstatic.get_bits(start, length);
      REQUIRE(random_bits == random_bits2);
    }
    bitsdynamic.set_bits(start, length, random_bits);
    chunk_t random_bits2 = bitsdynamic.get_bits(start, length);
    REQUIRE(random_bits == random_bits2);
  }
}

TEST_CASE("bitset", "[bits]") {
  Log.out("Testing bitset");
  
  test_bitset_get_bits<uint8_t>();
  test_bitset_get_bits<uint16_t>();
  test_bitset_get_bits<uint32_t>();
  test_bitset_get_bits<uint64_t>();

  test_bitset_set_bits<uint8_t>();
  test_bitset_set_bits<uint16_t>();
  test_bitset_set_bits<uint32_t>();
  test_bitset_set_bits<uint64_t>();
}
