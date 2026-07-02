// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <sstream>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/states/inner.hpp>

using namespace xdiag;

TEST_CASE("product_state", "[states]") try {

  // Test product state for Electron block
  for (int i = 0; i < 10; ++i) {
    for (int nsites = 1; nsites < 9; ++nsites) {

      std::vector<int> ps = {0, 1, 2, 3, 0, 1, 2, 3};
      std::vector<int> pss;
      std::sample(ps.begin(), ps.end(), std::back_inserter(pss), nsites,
                  std::mt19937(i));
      auto pstate = ProductState(pss);
      // XDIAG_SHOW(pstate);

      auto block = Electron(nsites);
      auto psi = State(block, true);
      fill(psi, pstate);

      int nup = 0;
      int ndn = 0;
      for (int i = 0; i < nsites; ++i) {
        int64_t p = pstate[i];
        auto sz = Op("Sz", i);
        auto szc = inner(sz, psi);

        auto n = Op("Ntot", i);
        auto nc = inner(n, psi);

        if (p == 0) {
          REQUIRE(isapprox(szc, 0.));
          REQUIRE(isapprox(nc, 0.));
        } else if (p == 1) {
          REQUIRE(isapprox(szc, 0.5));
          REQUIRE(isapprox(nc, 1.0));
          ++nup;
        } else if (p == 2) {
          REQUIRE(isapprox(szc, -0.5));
          REQUIRE(isapprox(nc, 1.0));
          ++ndn;
        } else if (p == 3) {
          REQUIRE(isapprox(szc, 0.0));
          REQUIRE(isapprox(nc, 2.0));
          ++nup;
          ++ndn;
        }
      }

      auto block2 = Electron(nsites, nup, ndn);
      auto psi2 = State(block2);
      fill(psi2, pstate);

      for (int i = 0; i < nsites; ++i) {
        int p = pstate[i];
        auto sz = Op("Sz", i);
        auto szc = inner(sz, psi2);

        auto n = Op("Ntot", i);
        auto nc = inner(n, psi2);

        if (p == 0) {
          REQUIRE(isapprox(szc, 0.));
          REQUIRE(isapprox(nc, 0.));
        } else if (p == 1) {
          REQUIRE(isapprox(szc, 0.5));
          REQUIRE(isapprox(nc, 1.0));
        } else if (p == 2) {
          REQUIRE(isapprox(szc, -0.5));
          REQUIRE(isapprox(nc, 1.0));
        } else if (p == 3) {
          REQUIRE(isapprox(szc, 0.0));
          REQUIRE(isapprox(nc, 2.0));
        }
      }
    }
  }

  // Test product state for tJ block
  for (int i = 0; i < 10; ++i) {
    for (int nsites = 1; nsites < 9; ++nsites) {

      std::vector<int> ps = {0, 1, 2, 0, 1, 2, 0, 1, 2};
      std::vector<int> pss;
      std::sample(ps.begin(), ps.end(), std::back_inserter(pss), nsites,
                  std::mt19937(i));
      auto pstate = ProductState(pss);
      // XDIAG_SHOW(pstate);

      int nup = 0;
      int ndn = 0;
      for (int i = 0; i < nsites; ++i) {
        int64_t p = pstate[i];

        if (p == 1) {
          ++nup;
        } else if (p == 2) {
          ++ndn;
        }
      }

      auto block2 = tJ(nsites, nup, ndn);
      auto psi2 = State(block2);
      fill(psi2, pstate);
      for (int i = 0; i < nsites; ++i) {
        int64_t p = pstate[i];
        auto sz = Op("Sz", i);
        auto szc = inner(sz, psi2);

        auto n = Op("Ntot", i);
        auto nc = inner(n, psi2);

        if (p == 0) {
          REQUIRE(isapprox(szc, 0.));
          REQUIRE(isapprox(nc, 0.));
        } else if (p == 1) {
          REQUIRE(isapprox(szc, 0.5));
          REQUIRE(isapprox(nc, 1.0));
        } else if (p == 2) {
          REQUIRE(isapprox(szc, -0.5));
          REQUIRE(isapprox(nc, 1.0));
        }
      }
    }
  }

  // Test product state for Spinhalf block
  for (int i = 0; i < 10; ++i) {
    for (int nsites = 1; nsites < 9; ++nsites) {

      std::vector<int64_t> ps = {1, 0, 1, 0, 1, 0, 1, 0};
      std::vector<int64_t> pss;
      std::sample(ps.begin(), ps.end(), std::back_inserter(pss), nsites,
                  std::mt19937(i));
      auto pstate = ProductState(pss);
      // XDIAG_SHOW(pstate);

      int nup = 0;
      auto block = Spinhalf(nsites);
      auto psi = State(block);
      fill(psi, pstate);
      for (int i = 0; i < nsites; ++i) {
        int64_t p = pstate[i];

        auto sz = Op("Sz", i);
        auto szc = inner(sz, psi);

        if (p == 1) {
          REQUIRE(isapprox(szc, 0.5));
          ++nup;
        } else if (p == 0) {
          REQUIRE(isapprox(szc, -0.5));
        }
      }

      auto block2 = Spinhalf(nsites, nup);
      auto psi2 = State(block2);
      fill(psi2, pstate);
      for (int i = 0; i < nsites; ++i) {
        int64_t p = pstate[i];
        auto sz = Op("Sz", i);
        auto szc = inner(sz, psi2);

        if (p == 1) {
          REQUIRE(isapprox(szc, 0.5));
        } else if (p == 0) {
          REQUIRE(isapprox(szc, -0.5));
        }
      }
    }
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
}

TEST_CASE("product_state_api", "[states]") try {
  using namespace xdiag;

  // --- Default construction is empty ---
  {
    ProductState p;
    REQUIRE(p.size() == 0);
    REQUIRE(p.nsites() == 0);
    REQUIRE(size(p) == 0);
    REQUIRE(nsites(p) == 0);
    REQUIRE(p.begin() == p.end());
  }

  // --- Sized construction is zero-initialized ---
  {
    ProductState p(3);
    REQUIRE(p.size() == 3);
    REQUIRE(p.nsites() == 3);
    for (int64_t i = 0; i < p.size(); ++i) {
      REQUIRE(p[i] == 0);
    }
  }

  // --- Construction from int64 and int32 vectors agree ---
  {
    std::vector<int64_t> v64 = {0, 1, 2, 1};
    std::vector<int32_t> v32 = {0, 1, 2, 1};
    ProductState p64(v64);
    ProductState p32(v32);
    REQUIRE(p64 == p32);
    REQUIRE(!(p64 != p32));
  }

  // --- Mutation via operator[] and push_back ---
  {
    ProductState p(2);
    p[0] = 1;
    p[1] = 2;
    p.push_back(3);
    REQUIRE(p.size() == 3);
    REQUIRE(p[0] == 1);
    REQUIRE(p[1] == 2);
    REQUIRE(p[2] == 3);

    // const operator[] through a const reference
    ProductState const &cp = p;
    REQUIRE(cp[2] == 3);
  }

  // --- Iteration yields the local states in order ---
  {
    std::vector<int64_t> v = {5, 6, 7};
    ProductState p(v);
    std::vector<int64_t> collected(p.begin(), p.end());
    REQUIRE(collected == v);
  }

  // --- Equality / inequality ---
  {
    ProductState a(std::vector<int64_t>{1, 0, 1});
    ProductState b(std::vector<int64_t>{1, 0, 1});
    ProductState c(std::vector<int64_t>{1, 1, 1});
    REQUIRE(a == b);
    REQUIRE(a != c);
    ProductState d(std::vector<int64_t>{1, 0}); // different length
    REQUIRE(a != d);
  }

  // --- Printing / to_string prints entries in reverse (MSB-first) order ---
  {
    ProductState p(std::vector<int64_t>{1, 2, 3});
    REQUIRE(to_string(p) == "3 2 1");

    std::ostringstream oss;
    oss << p;
    REQUIRE(oss.str() == "3 2 1");

    ProductState single(std::vector<int64_t>{7});
    REQUIRE(to_string(single) == "7");

    ProductState empty;
    REQUIRE(to_string(empty).empty());
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
}
