// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <iostream>
#include <vector>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/logger.hpp>

template <class bit_t>
void test_subsets(){
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n=0; n<8; ++n)
      {
	Subsets<bit_t> subs(n);
	REQUIRE(n == subs.n());

        int64_t ctr=0;
	bit_t current=0;
	for (auto sub : subs)
	  {
	    if (ctr != 0) REQUIRE(sub > current);
	    current = sub;
	    ++ctr;
	    REQUIRE(sub < ((bit_t)1 << n));
	  }
	REQUIRE(ctr == subs.size());
      }
}

template <class bit_t>
void test_subsets_random_access() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n = 0; n < 8; ++n) {
    Subsets<bit_t> subs(n);

    // operator[]: subs[i] == (bit_t)i
    for (int64_t i = 0; i < subs.size(); ++i)
      REQUIRE(subs[i] == (bit_t)i);

    // index: round-trip index(subs[i]) == i
    for (int64_t i = 0; i < subs.size(); ++i)
      REQUIRE(subs.index(subs[i]) == i);
  }
}

template <class bit_t>
void test_subsets_iterator_advance() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n = 1; n < 8; ++n) {
    Subsets<bit_t> subs(n);

    // Collect all elements via sequential iteration
    std::vector<bit_t> elems;
    for (auto s : subs)
      elems.push_back(s);

    // operator+: begin() + i should match elems[i]
    for (int64_t i = 0; i < subs.size(); ++i)
      REQUIRE(*(subs.begin() + i) == elems[i]);

    // operator+=: advance a single iterator step by step
    auto it = subs.begin();
    for (int64_t i = 0; i < subs.size() - 1; ++i) {
      it += 1;
      REQUIRE(*it == elems[i + 1]);
    }

    // operator+=: jump by larger steps
    if (subs.size() >= 4) {
      auto it2 = subs.begin();
      it2 += 3;
      REQUIRE(*it2 == elems[3]);
    }
  }
}

TEST_CASE( "subsets", "[combinatorics]" ) {

  SECTION("iteration") {
    xdiag::Log("Testing Subsets - iteration");
    test_subsets<uint32_t>();
    test_subsets<uint64_t>();
  }

  SECTION("random access and index") {
    xdiag::Log("Testing Subsets - random access and index");
    test_subsets_random_access<uint32_t>();
    test_subsets_random_access<uint64_t>();
  }

  SECTION("iterator advance (+ and +=)") {
    xdiag::Log("Testing Subsets - iterator advance");
    test_subsets_iterator_advance<uint32_t>();
    test_subsets_iterator_advance<uint64_t>();
  }
}
