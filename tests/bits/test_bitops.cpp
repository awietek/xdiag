// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <limits>
#include <random>

#include <xdiag/bits/ctz.hpp>
#include <xdiag/bits/extract_deposit.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/timing.hpp>

using namespace xdiag;
using namespace xdiag::bits;

template <typename bit_t> void test_bitops() {

  // constexpr int N = 10000000;
  // bool do_timing = true;

  constexpr int N = 100;
  bool do_timing = false;

  std::mt19937_64 gen;

  std::vector<bit_t> v1(N);
  std::vector<bit_t> v2(N);
  for (int i = 0; i < N; ++i) {
    v1[i] = gen();
    v2[i] = gen();
  }

  // Compare popcounts
  {
    std::vector<int> popcounts1(N);
    std::vector<int> popcounts2(N);

    auto time = rightnow();
    for (int i = 0; i < N; ++i) {
      popcounts1[i] = swar_popcount(v1[i]);
    }
    if (do_timing)
      timing(time, rightnow(), "popcount (slow)", 1);

    time = rightnow();
    for (int i = 0; i < N; ++i) {
      popcounts2[i] = popcount(v1[i]);
    }
    if (do_timing)
      timing(time, rightnow(), "popcount (fast)", 1);

    for (int i = 0; i < N; ++i) {
      REQUIRE(popcounts1[i] == popcounts2[i]);
    }
  }

  // Compare extract
  {
    std::vector<int> pexts1(N);
    std::vector<int> pexts2(N);

    auto time = rightnow();
    for (int i = 0; i < N; ++i) {
      pexts1[i] = pext_fallback(v1[i], v2[i]);
    }
    if (do_timing)
      timing(time, rightnow(), "pext (slow)", 1);

    time = rightnow();
    for (int i = 0; i < N; ++i) {
      pexts2[i] = extract(v1[i], v2[i]);
    }
    if (do_timing)
      timing(time, rightnow(), "pext (fast)", 1);

    for (int i = 0; i < N; ++i) {
      REQUIRE(pexts1[i] == pexts2[i]);
    }
  }

  // Compare deposit
  {
    std::vector<int> pdeps1(N);
    std::vector<int> pdeps2(N);

    auto time = rightnow();
    for (int i = 0; i < N; ++i) {
      pdeps1[i] = pdep_fallback(v1[i], v2[i]);
    }
    if (do_timing)
      timing(time, rightnow(), "pdep (slow)", 1);

    time = rightnow();
    for (int i = 0; i < N; ++i) {
      pdeps2[i] = deposit(v1[i], v2[i]);
    }
    if (do_timing)
      timing(time, rightnow(), "pdep (fast)", 1);

    for (int i = 0; i < N; ++i) {
      REQUIRE(pdeps1[i] == pdeps2[i]);
    }
  }

  // Compare gbit
  {
    int n = 7;

    std::vector<int> gbits1(N);
    std::vector<int> gbits2(N);

    auto time = rightnow();
    for (int i = 0; i < N; ++i) {
      gbits1[i] = (v1[i] >> n) & 1;
    }
    if (do_timing)
      timing(time, rightnow(), "gbit (slow)", 1);

    time = rightnow();
    for (int i = 0; i < N; ++i) {
      gbits2[i] = get(v1[i], n);
    }
    if (do_timing)
      timing(time, rightnow(), "gbit (fast)", 1);

    for (int i = 0; i < N; ++i) {
      REQUIRE(gbits1[i] == gbits2[i]);
    }
  }

  // Compare gbits
  {
    int n = 13;
    int l = 17;

    std::vector<int> gbitss1(N);
    std::vector<int> gbitss2(N);

    auto time = rightnow();
    for (int i = 0; i < N; ++i) {
      gbitss1[i] = (v1[i] >> n) & ((1 << l) - 1);
    }
    if (do_timing)
      timing(time, rightnow(), "gbitss (slow)", 1);

    time = rightnow();
    for (int i = 0; i < N; ++i) {
      gbitss2[i] = get_range(v1[i], n, l);
    }
    if (do_timing)
      timing(time, rightnow(), "gbitss (fast)", 1);

    for (int i = 0; i < N; ++i) {
      REQUIRE(gbitss1[i] == gbitss2[i]);
    }
  }
}

TEST_CASE("bitops", "[bitops]") {
  Log.set_verbosity(1);
  Log("Testing bitops<uint32_t>");
  test_bitops<uint32_t>();
  Log("Testing bitops<uint64_t>");
  test_bitops<uint64_t>();
  Log.set_verbosity(0);
}

// trailing zeros by plain bit scanning (x != 0)
template <typename bit_t> static int ctz_reference(bit_t x) {
  int n = 0;
  while (!((x >> n) & 1)) {
    ++n;
  }
  return n;
}

// ctz (builtin, or the fallback on other compilers) and the SWAR fallback,
// which is tested explicitly as gcc / clang never use it
template <typename bit_t> void test_ctz(int (*swar_ctz)(bit_t)) {
  constexpr int width = std::numeric_limits<bit_t>::digits;

  // bit i alone, and bit i with all bits above set: the higher bits must not
  // change the result
  for (int i = 0; i < width; ++i) {
    bit_t x = bit_t(1) << i;
    bit_t y = ~bit_t(0) << i;
    REQUIRE(ctz(x) == i);
    REQUIRE(ctz(y) == i);
    REQUIRE(swar_ctz(x) == i);
    REQUIRE(swar_ctz(y) == i);
  }

  // random values with at least k trailing zeros, k running over all widths
  std::mt19937_64 gen(42);
  for (int k = 0; k < 1000; ++k) {
    bit_t x = (bit_t)gen() << (k % width);
    if (x == 0) {
      continue;
    }
    REQUIRE(ctz(x) == ctz_reference(x));
    REQUIRE(swar_ctz(x) == ctz_reference(x));
  }

  // the fallback is defined for zero and gives the full width
  REQUIRE(swar_ctz(bit_t(0)) == width);
}

TEST_CASE("ctz", "[bitops]") {
  Log("Testing ctz<uint32_t>");
  test_ctz<uint32_t>(swar_ctz_32);
  Log("Testing ctz<uint64_t>");
  test_ctz<uint64_t>(swar_ctz_64);
}
