#include "../catch.hpp"

#include <hydra/all.h>

template <class bit_t> void test_state_spinhalf() {
  using namespace hydra;

  bit_t max_bits = 16;
  for (bit_t b1 = 0; b1 < max_bits; ++b1) {
    auto s1 = state_spinhalf<bit_t>({b1});

    for (bit_t b2 = 0; b2 < max_bits; ++b2) {
      auto s2 = state_spinhalf<bit_t>({b2});

      if (b1 == b2)
        REQUIRE(s1 == s2);
      else
        REQUIRE(s1 != s2);

      if (b1 < b2)
        REQUIRE(s1 < s2);
      else
        REQUIRE(s1 >= s2);

      if (b1 > b2)
        REQUIRE(s1 > s2);
      else
        REQUIRE(s1 <= s2);
    }
    REQUIRE(QN(s1).n_up == utils::popcnt(b1));

    for (int shift=0; shift < 4; ++shift)
      {
	auto s1s = s1 << shift;
	REQUIRE(s1s.spins == s1.spins << shift);
	REQUIRE(s1 == s1s >> shift);

	s1s = s1;
	s1s <<= shift;
	REQUIRE(s1s.spins == s1.spins << shift);
	s1s >>= shift;
	REQUIRE(s1s == s1);
      }

    int min_bit = 1;
    int max_bit = 4;
    auto ssub = sitevals(s1, max_bit - min_bit, min_bit);
    
    auto snew1 = state_spinhalf<bit_t>({0});
    auto snew2 = state_spinhalf<bit_t>({0});
    for (int bit=min_bit; bit<max_bit; ++bit)
      {
	auto sv = siteval(s1, bit);
	snew1 = snew1 | (sv << (bit - min_bit));
	snew2 |= sv << (bit - min_bit);
      }
    REQUIRE(ssub == snew1);
    REQUIRE(ssub == snew2);    
  }
}

TEST_CASE("states", "[states]") {
  test_state_spinhalf<uint16>();
  test_state_spinhalf<uint32>();
  test_state_spinhalf<uint64>();
}
