// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/symmetries/action/site_permutation.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/utils/logger.hpp>

TEST_CASE("site_permutation", "[symmetries]") try {
  using namespace xdiag::symmetries;
  using namespace xdiag;
  using namespace xdiag::bits;
  Log("Test SitePermutation");

  // Basic accessors
  for (int64_t n = 2; n < 7; ++n) {
    auto sp = SitePermutation(cyclic_group(n));
    REQUIRE(sp.size() == n);
    REQUIRE(sp.nsites() == n);
  }

  // operator== / !=
  {
    auto sp1 = SitePermutation(cyclic_group(4));
    auto sp2 = SitePermutation(cyclic_group(4));
    auto sp3 = SitePermutation(cyclic_group(5));
    REQUIRE(sp1 == sp2);
    REQUIRE(sp1 != sp3);
  }

  // apply<uint16_t>: cyclic group on 4 sites, sym k shifts bit i -> site (i+k)%4
  {
    auto sp = SitePermutation(cyclic_group(4));
    REQUIRE(sp.apply<uint16_t>(0, 0b0001) == 0b0001); // identity
    REQUIRE(sp.apply<uint16_t>(1, 0b0001) == 0b0010); // bit 0 -> site 1
    REQUIRE(sp.apply<uint16_t>(2, 0b0001) == 0b0100); // bit 0 -> site 2
    REQUIRE(sp.apply<uint16_t>(3, 0b0001) == 0b1000); // bit 0 -> site 3
    REQUIRE(sp.apply<uint16_t>(1, 0b0101) == 0b1010); // bits {0,2} -> {1,3}
  }

  // apply(inv(sym), apply(sym, bits)) == bits for all states and symmetries
  {
    auto group = cyclic_group(4);
    auto sp = SitePermutation(group);
    for (uint16_t state = 0; state < 16; ++state) {
      for (int64_t sym = 0; sym < group.size(); ++sym) {
        auto permuted = sp.apply<uint16_t>(sym, state);
        auto recovered = sp.apply<uint16_t>(group.inv(sym), permuted);
        REQUIRE(recovered == state);
      }
    }
  }

  // Homomorphism: apply(s2, apply(s1, b)) == apply(multiply(s1,s2), b)
  {
    auto group = cyclic_group(4);
    auto sp = SitePermutation(group);
    for (uint16_t state = 0; state < 16; ++state) {
      for (int64_t s1 = 0; s1 < group.size(); ++s1) {
        for (int64_t s2 = 0; s2 < group.size(); ++s2) {
          auto composed = sp.apply<uint16_t>(s2, sp.apply<uint16_t>(s1, state));
          auto direct = sp.apply<uint16_t>(group.multiply(s1, s2), state);
          REQUIRE(composed == direct);
        }
      }
    }
  }

  // apply with BitArray<uint64_t, 1>: 1-bit-per-site storage
  {
    auto group = cyclic_group(4);
    auto sp = SitePermutation(group);
    BitArray<uint64_t, 1> bits_in;
    bits_in.set(0, 1); // only bit 0 set
    // sym k should move the set bit from site 0 to site k
    for (int64_t sym = 0; sym < group.size(); ++sym) {
      auto bits_out = sp.apply(sym, bits_in);
      for (int64_t j = 0; j < 4; ++j)
        REQUIRE(bits_out.get(j) == (j == sym ? 1 : 0));
    }
  }

  // apply with BitArray<uint64_t, 2>: 2-bit values per site
  {
    auto group = cyclic_group(4);
    auto sp = SitePermutation(group);
    // assign distinct 2-bit values to each site: 3, 2, 1, 0
    BitArray<uint64_t, 2> bits_in;
    for (int64_t i = 0; i < 4; ++i)
      bits_in.set(i, 3 - i); // site 0=3, site 1=2, site 2=1, site 3=0
    // sym k: value at site i moves to site (i+k)%4
    for (int64_t sym = 0; sym < group.size(); ++sym) {
      auto bits_out = sp.apply(sym, bits_in);
      for (int64_t i = 0; i < 4; ++i)
        REQUIRE(bits_out.get((i + sym) % 4) == bits_in.get(i));
    }
  }

  // apply with BitArray<uint64_t, 3>: 3-bit values per site
  {
    auto group = cyclic_group(4);
    auto sp = SitePermutation(group);
    // assign distinct 3-bit values to each site: 5, 3, 7, 1
    BitArray<uint64_t, 3> bits_in;
    std::array<int64_t, 4> vals = {5, 3, 7, 1};
    for (int64_t i = 0; i < 4; ++i)
      bits_in.set(i, vals[i]);
    for (int64_t sym = 0; sym < group.size(); ++sym) {
      auto bits_out = sp.apply(sym, bits_in);
      for (int64_t i = 0; i < 4; ++i)
        REQUIRE(bits_out.get((i + sym) % 4) == bits_in.get(i));
    }
  }

  // apply with BitsetDynamic: must be constructed with explicit size
  {
    auto group = cyclic_group(4);
    auto sp = SitePermutation(group);
    BitsetDynamic bits_in(4); // 4-bit dynamic bitset
    bits_in.set(0);           // set site 0
    REQUIRE(sp.apply(0, bits_in).test(0) == true);
    REQUIRE(sp.apply(1, bits_in).test(1) == true);
    REQUIRE(sp.apply(2, bits_in).test(2) == true);
    REQUIRE(sp.apply(3, bits_in).test(3) == true);
    // inverse recovery
    for (int64_t sym = 0; sym < group.size(); ++sym) {
      auto permuted = sp.apply(sym, bits_in);
      auto recovered = sp.apply(group.inv(sym), permuted);
      REQUIRE(recovered == bits_in);
    }
  }

  // apply with BitsetStatic1 (Bitset<uint64_t, 1>, 64-bit storage): bit permutation
  {
    auto group = cyclic_group(4);
    auto sp = SitePermutation(group);
    // single bit at site 0: sym k moves it to site k
    auto bits_in = make_bitset<uint64_t, 1>(0b0001);
    REQUIRE(to_uint64(sp.apply(0, bits_in)) == 0b0001);
    REQUIRE(to_uint64(sp.apply(1, bits_in)) == 0b0010);
    REQUIRE(to_uint64(sp.apply(2, bits_in)) == 0b0100);
    REQUIRE(to_uint64(sp.apply(3, bits_in)) == 0b1000);
    // multi-bit: {0,2} -> {1,3} under sym 1
    REQUIRE(to_uint64(sp.apply(1, make_bitset<uint64_t, 1>(0b0101))) == 0b1010);
  }

  // apply with BitsetStatic2 (Bitset<uint64_t, 2>, 128-bit storage): bit permutation
  {
    auto group = cyclic_group(4);
    auto sp = SitePermutation(group);
    auto bits_in = make_bitset<uint64_t, 2>(0b0001);
    REQUIRE(to_uint64(sp.apply(0, bits_in)) == 0b0001);
    REQUIRE(to_uint64(sp.apply(1, bits_in)) == 0b0010);
    REQUIRE(to_uint64(sp.apply(2, bits_in)) == 0b0100);
    REQUIRE(to_uint64(sp.apply(3, bits_in)) == 0b1000);
  }

  Log("done");
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
  throw;
}
