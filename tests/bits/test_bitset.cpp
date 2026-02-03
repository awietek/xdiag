// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/bits/bitset.hpp>
#include <xdiag/utils/logger.hpp>

#include <limits>
#include <random>

using namespace xdiag;
using namespace xdiag::bits;

// Helper function to convert uint64_t to bitset
template <typename chunk_t, int64_t nchunks>
static Bitset<chunk_t, nchunks> make_bitset(uint64_t value) {
  Bitset<chunk_t, nchunks> bits(64);
  for (int i = 0; i < 64; ++i) {
    if (value & (1ULL << i)) {
      bits.set(i);
    }
  }
  return bits;
}

// Helper function to convert bitset to uint64_t
template <typename chunk_t, int64_t nchunks>
static uint64_t to_uint64(Bitset<chunk_t, nchunks> const &bits) {
  uint64_t result = 0;
  for (int i = 0; i < 64; ++i) {
    if (bits.test(i)) {
      result |= (1ULL << i);
    }
  }
  return result;
}

template <typename chunk_t> void test_bitset_constructor() {
  constexpr int chunkdigits = std::numeric_limits<chunk_t>::digits;
  constexpr int nchunks = 64 / chunkdigits;

  // Test dynamic bitset
  {
    auto bits = Bitset<chunk_t, 0>(64);
    for (int i = 0; i < 64; ++i) {
      REQUIRE(bits.test(i) == false);
    }
  }

  // Test static bitset
  if constexpr (chunkdigits > 8) {
    auto bits = Bitset<chunk_t, nchunks>(64);
    for (int i = 0; i < 64; ++i) {
      REQUIRE(bits.test(i) == false);
    }
  }
}

template <typename chunk_t> void test_bitset_test_set_reset_flip() {
  constexpr int chunkdigits = std::numeric_limits<chunk_t>::digits;
  constexpr int nchunks = 64 / chunkdigits;

  std::mt19937 rng(42);
  std::uniform_int_distribution<int> dist(0, 63);

  for (int trial = 0; trial < 100; ++trial) {
    uint64_t value = 0;
    auto bits_dynamic = Bitset<chunk_t, 0>(64);
    Bitset<chunk_t, nchunks> bits_static(64);

    // Test set
    for (int i = 0; i < 20; ++i) {
      int pos = dist(rng);
      value |= (1ULL << pos);
      bits_dynamic.set(pos);
      if constexpr (chunkdigits > 8) {
        bits_static.set(pos);
      }
    }

    REQUIRE(to_uint64(bits_dynamic) == value);
    if constexpr (chunkdigits > 8) {
      REQUIRE(to_uint64(bits_static) == value);
    }

    // Test reset
    for (int i = 0; i < 10; ++i) {
      int pos = dist(rng);
      value &= ~(1ULL << pos);
      bits_dynamic.reset(pos);
      if constexpr (chunkdigits > 8) {
        bits_static.reset(pos);
      }
    }

    REQUIRE(to_uint64(bits_dynamic) == value);
    if constexpr (chunkdigits > 8) {
      REQUIRE(to_uint64(bits_static) == value);
    }

    // Test flip
    for (int i = 0; i < 10; ++i) {
      int pos = dist(rng);
      value ^= (1ULL << pos);
      bits_dynamic.flip(pos);
      if constexpr (chunkdigits > 8) {
        bits_static.flip(pos);
      }
    }

    REQUIRE(to_uint64(bits_dynamic) == value);
    if constexpr (chunkdigits > 8) {
      REQUIRE(to_uint64(bits_static) == value);
    }

    // Test set with value
    for (int i = 0; i < 10; ++i) {
      int pos = dist(rng);
      bool val = (i % 2 == 0);
      if (val) {
        value |= (1ULL << pos);
      } else {
        value &= ~(1ULL << pos);
      }
      bits_dynamic.set(pos, val);
      if constexpr (chunkdigits > 8) {
        bits_static.set(pos, val);
      }
    }

    REQUIRE(to_uint64(bits_dynamic) == value);
    if constexpr (chunkdigits > 8) {
      REQUIRE(to_uint64(bits_static) == value);
    }
  }
}

template <typename chunk_t> void test_bitset_bitwise_ops() {
  constexpr int chunkdigits = std::numeric_limits<chunk_t>::digits;
  constexpr int nchunks = 64 / chunkdigits;

  std::mt19937 rng(42);
  std::uniform_int_distribution<uint64_t> dist(
      std::numeric_limits<uint64_t>::min(),
      std::numeric_limits<uint64_t>::max());

  for (int trial = 0; trial < 100; ++trial) {
    uint64_t a = dist(rng);
    uint64_t b = dist(rng);

    auto bits_a = make_bitset<chunk_t, 0>(a);
    auto bits_b = make_bitset<chunk_t, 0>(b);

    // Test AND
    auto bits_and = bits_a & bits_b;
    REQUIRE(to_uint64(bits_and) == (a & b));

    // Test OR
    auto bits_or = bits_a | bits_b;
    REQUIRE(to_uint64(bits_or) == (a | b));

    // Test XOR
    auto bits_xor = bits_a ^ bits_b;
    REQUIRE(to_uint64(bits_xor) == (a ^ b));

    // Test NOT
    auto bits_not = ~bits_a;
    REQUIRE(to_uint64(bits_not) == ~a);

    // Test compound assignments
    auto bits_temp = bits_a;
    bits_temp &= bits_b;
    REQUIRE(to_uint64(bits_temp) == (a & b));

    bits_temp = bits_a;
    bits_temp |= bits_b;
    REQUIRE(to_uint64(bits_temp) == (a | b));

    bits_temp = bits_a;
    bits_temp ^= bits_b;
    REQUIRE(to_uint64(bits_temp) == (a ^ b));

    // Test static bitset
    if constexpr (chunkdigits > 8) {
      auto bits_a_s = make_bitset<chunk_t, nchunks>(a);
      auto bits_b_s = make_bitset<chunk_t, nchunks>(b);

      REQUIRE(to_uint64(bits_a_s & bits_b_s) == (a & b));
      REQUIRE(to_uint64(bits_a_s | bits_b_s) == (a | b));
      REQUIRE(to_uint64(bits_a_s ^ bits_b_s) == (a ^ b));
      REQUIRE(to_uint64(~bits_a_s) == ~a);
    }
  }
}

template <typename chunk_t> void test_bitset_shift_ops() {
  constexpr int chunkdigits = std::numeric_limits<chunk_t>::digits;
  constexpr int nchunks = 64 / chunkdigits;

  std::mt19937 rng(42);
  std::uniform_int_distribution<uint64_t> dist(
      std::numeric_limits<uint64_t>::min(),
      std::numeric_limits<uint64_t>::max());
  std::uniform_int_distribution<int> shift_dist(0, 63);

  for (int trial = 0; trial < 100; ++trial) {
    uint64_t value = dist(rng);
    int shift = shift_dist(rng);

    auto bits = make_bitset<chunk_t, 0>(value);

    // Test left shift
    auto bits_lshift = bits << shift;
    REQUIRE(to_uint64(bits_lshift) == (value << shift));

    // Test right shift
    auto bits_rshift = bits >> shift;
    REQUIRE(to_uint64(bits_rshift) == (value >> shift));

    // Test compound left shift
    auto bits_temp = bits;
    bits_temp <<= shift;
    REQUIRE(to_uint64(bits_temp) == (value << shift));

    // Test compound right shift
    bits_temp = bits;
    bits_temp >>= shift;
    REQUIRE(to_uint64(bits_temp) == (value >> shift));

    // Test static bitset
    if constexpr (chunkdigits > 8) {
      auto bits_s = make_bitset<chunk_t, nchunks>(value);
      REQUIRE(to_uint64(bits_s << shift) == (value << shift));
      REQUIRE(to_uint64(bits_s >> shift) == (value >> shift));
    }
  }
}

template <typename chunk_t> void test_bitset_predicates() {
  constexpr int chunkdigits = std::numeric_limits<chunk_t>::digits;
  constexpr int nchunks = 64 / chunkdigits;

  // Test all()
  {
    auto bits = Bitset<chunk_t, 0>(64);
    REQUIRE(bits.none() == true);
    REQUIRE(bits.any() == false);
    REQUIRE(bits.all() == false);
    REQUIRE(bits.count() == 0);

    for (int i = 0; i < 64; ++i) {
      bits.set(i);
    }
    REQUIRE(bits.all() == true);
    REQUIRE(bits.any() == true);
    REQUIRE(bits.none() == false);
    REQUIRE(bits.count() == 64);
  }

  // Test any() and none()
  {
    auto bits = Bitset<chunk_t, 0>(64);
    REQUIRE(bits.none() == true);
    REQUIRE(bits.any() == false);

    bits.set(32);
    REQUIRE(bits.none() == false);
    REQUIRE(bits.any() == true);
  }

  // Test count() with random values
  std::mt19937 rng(42);
  std::uniform_int_distribution<uint64_t> dist(
      std::numeric_limits<uint64_t>::min(),
      std::numeric_limits<uint64_t>::max());

  for (int trial = 0; trial < 100; ++trial) {
    uint64_t value = dist(rng);
    auto bits = make_bitset<chunk_t, 0>(value);

    int expected_count = 0;
    for (int i = 0; i < 64; ++i) {
      if (value & (1ULL << i)) {
        expected_count++;
      }
    }

    REQUIRE(bits.count() == expected_count);
    REQUIRE(bits.any() == (expected_count > 0));
    REQUIRE(bits.none() == (expected_count == 0));
    REQUIRE(bits.all() == (expected_count == 64));

    // Test static bitset
    if constexpr (chunkdigits > 8) {
      auto bits_s = make_bitset<chunk_t, nchunks>(value);
      REQUIRE(bits_s.count() == expected_count);
      REQUIRE(bits_s.any() == (expected_count > 0));
      REQUIRE(bits_s.none() == (expected_count == 0));
      REQUIRE(bits_s.all() == (expected_count == 64));
    }
  }
}

template <typename chunk_t> void test_bitset_comparison() {
  constexpr int chunkdigits = std::numeric_limits<chunk_t>::digits;
  constexpr int nchunks = 64 / chunkdigits;

  std::mt19937 rng(42);
  std::uniform_int_distribution<uint64_t> dist(
      std::numeric_limits<uint64_t>::min(),
      std::numeric_limits<uint64_t>::max());

  for (int trial = 0; trial < 100; ++trial) {
    uint64_t a = dist(rng);
    uint64_t b = dist(rng);

    auto bits_a = make_bitset<chunk_t, 0>(a);
    auto bits_b = make_bitset<chunk_t, 0>(b);
    auto bits_a_copy = make_bitset<chunk_t, 0>(a);

    REQUIRE((bits_a == bits_b) == (a == b));
    REQUIRE((bits_a != bits_b) == (a != b));
    REQUIRE(bits_a == bits_a_copy);
    REQUIRE(!(bits_a != bits_a_copy));

    // Test static bitset
    if constexpr (chunkdigits > 8) {
      auto bits_a_s = make_bitset<chunk_t, nchunks>(a);
      auto bits_b_s = make_bitset<chunk_t, nchunks>(b);
      auto bits_a_s_copy = make_bitset<chunk_t, nchunks>(a);

      REQUIRE((bits_a_s == bits_b_s) == (a == b));
      REQUIRE((bits_a_s != bits_b_s) == (a != b));
      REQUIRE(bits_a_s == bits_a_s_copy);
      REQUIRE(!(bits_a_s != bits_a_s_copy));
    }
  }
}

template <typename chunk_t> void test_bitset_get_range() {
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
      auto bits = make_bitset<chunk_t, nchunks>(random_uint64);
      chunk_t bits_filtered = bits.get_range(start, length);
      REQUIRE(random_uint64_filtered == bits_filtered);
    }
    auto bits = make_bitset<chunk_t, 0>(random_uint64);
    chunk_t bits_filtered = bits.get_range(start, length);
    REQUIRE(random_uint64_filtered == bits_filtered);
  }
}

template <typename chunk_t> void test_bitset_set_range() {
  constexpr int chunkdigits = std::numeric_limits<chunk_t>::digits;
  constexpr int nchunks = 64 / chunkdigits;
  auto bitsstatic = Bitset<chunk_t, nchunks>();
  auto bitsdynamic = Bitset<chunk_t>(64);

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
      bitsstatic.set_range(start, length, random_bits);
      chunk_t random_bits2 = bitsstatic.get_range(start, length);
      REQUIRE(random_bits == random_bits2);
    }
    bitsdynamic.set_range(start, length, random_bits);
    chunk_t random_bits2 = bitsdynamic.get_range(start, length);
    REQUIRE(random_bits == random_bits2);
  }
}

TEST_CASE("bitset", "[bits]") {
  Log.out("Testing bitset");

  SECTION("constructor") {
    test_bitset_constructor<uint8_t>();
    test_bitset_constructor<uint16_t>();
    test_bitset_constructor<uint32_t>();
    test_bitset_constructor<uint64_t>();
  }

  SECTION("test_set_reset_flip") {
    test_bitset_test_set_reset_flip<uint8_t>();
    test_bitset_test_set_reset_flip<uint16_t>();
    test_bitset_test_set_reset_flip<uint32_t>();
    test_bitset_test_set_reset_flip<uint64_t>();
  }

  SECTION("get_range") {
    test_bitset_get_range<uint8_t>();
    test_bitset_get_range<uint16_t>();
    test_bitset_get_range<uint32_t>();
    test_bitset_get_range<uint64_t>();
  }

  SECTION("set_range") {
    test_bitset_set_range<uint8_t>();
    test_bitset_set_range<uint16_t>();
    test_bitset_set_range<uint32_t>();
    test_bitset_set_range<uint64_t>();
  }

  SECTION("bitwise_ops") {
    test_bitset_bitwise_ops<uint8_t>();
    test_bitset_bitwise_ops<uint16_t>();
    test_bitset_bitwise_ops<uint32_t>();
    test_bitset_bitwise_ops<uint64_t>();
  }

  SECTION("shift_ops") {
    test_bitset_shift_ops<uint8_t>();
    test_bitset_shift_ops<uint16_t>();
    test_bitset_shift_ops<uint32_t>();
    test_bitset_shift_ops<uint64_t>();
  }

  SECTION("predicates") {
    test_bitset_predicates<uint8_t>();
    test_bitset_predicates<uint16_t>();
    test_bitset_predicates<uint32_t>();
    test_bitset_predicates<uint64_t>();
  }

  SECTION("comparison") {
    test_bitset_comparison<uint8_t>();
    test_bitset_comparison<uint16_t>();
    test_bitset_comparison<uint32_t>();
    test_bitset_comparison<uint64_t>();
  }
}
