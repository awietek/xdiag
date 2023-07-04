#include "../catch.hpp"

#include <hydra/bitops/bitops.h>
#include <hydra/utils/timing.h>
#include <hydra/utils/logger.h>

#include <random>

using namespace hydra;
using namespace hydra::bitops;

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

  // Compare popcnts
  {
    std::vector<int> popcnts1(N);
    std::vector<int> popcnts2(N);

    auto time = rightnow();
    for (int i = 0; i < N; ++i) {
      popcnts1[i] = swar_popcnt(v1[i]);
    }
    if (do_timing)
      timing(time, rightnow(), "popcnt (slow)", 1);

    time = rightnow();
    for (int i = 0; i < N; ++i) {
      popcnts2[i] = popcnt(v1[i]);
    }
    if (do_timing)
      timing(time, rightnow(), "popcnt (fast)", 1);

    for (int i = 0; i < N; ++i) {
      REQUIRE(popcnts1[i] == popcnts2[i]);
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
      gbits2[i] = gbit(v1[i], n);
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
      gbitss2[i] = gbits(v1[i], n, l);
    }
    if (do_timing)
      timing(time, rightnow(), "gbitss (fast)", 1);

    for (int i = 0; i < N; ++i) {
      REQUIRE(gbitss1[i] == gbitss2[i]);
    }
  }
}

TEST_CASE("bitops", "[bitops]") {
  Log.out("Testing bitops");
  Log.set_verbosity(1);
  test_bitops<uint16_t>();
  test_bitops<uint32_t>();
  test_bitops<uint64_t>();
  Log.set_verbosity(0);
}
