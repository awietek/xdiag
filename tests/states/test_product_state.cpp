// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/algebra/isapprox.hpp>

using namespace xdiag;

TEST_CASE("product_state", "[states]") try {

  // Test product state for Electron block
  for (int i = 0; i < 10; ++i) {
    for (int nsites = 1; nsites < 9; ++nsites) {

      std::vector<std::string> ps = {"Emp", "Up", "Dn", "UpDn",
                                     "Emp", "Up", "Dn", "UpDn"};
      std::vector<std::string> pss;
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
        std::string p = pstate[i];
        auto sz = Op("Sz", i);
        auto szc = inner(sz, psi);

        auto n = Op("Ntot", i);
        auto nc = inner(n, psi);

        if (p == "Emp") {
          REQUIRE(isapprox(szc, 0.));
          REQUIRE(isapprox(nc, 0.));
        } else if (p == "Up") {
          REQUIRE(isapprox(szc, 0.5));
          REQUIRE(isapprox(nc, 1.0));
          ++nup;
        } else if (p == "Dn") {
          REQUIRE(isapprox(szc, -0.5));
          REQUIRE(isapprox(nc, 1.0));
          ++ndn;
        } else if (p == "UpDn") {
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
        std::string p = pstate[i];
        auto sz = Op("Sz", i);
        auto szc = inner(sz, psi2);

        auto n = Op("Ntot", i);
        auto nc = inner(n, psi2);

        if (p == "Emp") {
          REQUIRE(isapprox(szc, 0.));
          REQUIRE(isapprox(nc, 0.));
        } else if (p == "Up") {
          REQUIRE(isapprox(szc, 0.5));
          REQUIRE(isapprox(nc, 1.0));
        } else if (p == "Dn") {
          REQUIRE(isapprox(szc, -0.5));
          REQUIRE(isapprox(nc, 1.0));
        } else if (p == "UpDn") {
          REQUIRE(isapprox(szc, 0.0));
          REQUIRE(isapprox(nc, 2.0));
        }
      }
    }
  }

  // Test product state for tJ block
  for (int i = 0; i < 10; ++i) {
    for (int nsites = 1; nsites < 9; ++nsites) {

      std::vector<std::string> ps = {
          "Emp", "Up", "Dn", "Emp", "Up", "Dn", "Emp", "Up", "Dn",
      };
      std::vector<std::string> pss;
      std::sample(ps.begin(), ps.end(), std::back_inserter(pss), nsites,
                  std::mt19937(i));
      auto pstate = ProductState(pss);
      // XDIAG_SHOW(pstate);

      int nup = 0;
      int ndn = 0;
      for (int i = 0; i < nsites; ++i) {
        std::string p = pstate[i];

        if (p == "Up") {
          ++nup;
        } else if (p == "Dn") {
          ++ndn;
        }
      }

      auto block2 = tJ(nsites, nup, ndn);
      auto psi2 = State(block2);
      fill(psi2, pstate);
      for (int i = 0; i < nsites; ++i) {
        std::string p = pstate[i];
        auto sz = Op("Sz", i);
        auto szc = inner(sz, psi2);

        auto n = Op("Ntot", i);
        auto nc = inner(n, psi2);

        if (p == "Emp") {
          REQUIRE(isapprox(szc, 0.));
          REQUIRE(isapprox(nc, 0.));
        } else if (p == "Up") {
          REQUIRE(isapprox(szc, 0.5));
          REQUIRE(isapprox(nc, 1.0));
        } else if (p == "Dn") {
          REQUIRE(isapprox(szc, -0.5));
          REQUIRE(isapprox(nc, 1.0));
        }
      }
    }
  }

  // Test product state for Spinhalf block
  for (int i = 0; i < 10; ++i) {
    for (int nsites = 1; nsites < 9; ++nsites) {

      std::vector<std::string> ps = {"Up", "Dn", "Up", "Dn",
                                     "Up", "Dn", "Up", "Dn"};
      std::vector<std::string> pss;
      std::sample(ps.begin(), ps.end(), std::back_inserter(pss), nsites,
                  std::mt19937(i));
      auto pstate = ProductState(pss);
      // XDIAG_SHOW(pstate);

      int nup = 0;
      auto block = Spinhalf(nsites);
      auto psi = State(block);
      fill(psi, pstate);
      for (int i = 0; i < nsites; ++i) {
        std::string p = pstate[i];

        auto sz = Op("Sz", i);
        auto szc = inner(sz, psi);

        if (p == "Up") {
          REQUIRE(isapprox(szc, 0.5));
          ++nup;
        } else if (p == "Dn") {
          REQUIRE(isapprox(szc, -0.5));
        }
      }

      auto block2 = Spinhalf(nsites, nup);
      auto psi2 = State(block2);
      fill(psi2, pstate);
      for (int i = 0; i < nsites; ++i) {
        std::string p = pstate[i];
        auto sz = Op("Sz", i);
        auto szc = inner(sz, psi2);

        if (p == "Up") {
          REQUIRE(isapprox(szc, 0.5));
        } else if (p == "Dn") {
          REQUIRE(isapprox(szc, -0.5));
        }
      }
    }
  }

  // {
  //   for

  //     auto pstate = ProductState({"Up", "Dn", "UpDn", "Emp"});
  //   auto block = Electron(4, 2, 2);
  //   auto psi = StateReal(block, pstate);

  //   auto sz0 = Op("Sz", 0);
  //   complex sz0c = inner(sz0, psi);
  //   REQUIRE(isapprox(sz0c, 0.5));

  //   auto sz1 = Op("Sz", 1);
  //   complex sz1c = inner(sz1, psi);
  //   REQUIRE(isapprox(sz1c, -0.5));

  //   auto sz2 = Op("Sz", 2);
  //   complex sz2c = inner(sz2, psi);
  //   REQUIRE(isapprox(sz2c, 0.));

  //   auto sz3 = Op("Sz", 3);
  //   complex sz3c = inner(sz3, psi);
  //   REQUIRE(isapprox(sz3c, 0.));

  //   auto n0 = Op("Ntot", 0);
  //   complex n0c = inner(n0, psi);
  //   REQUIRE(isapprox(n0c, 1.0));

  //   auto n1 = Op("Ntot", 1);
  //   complex n1c = inner(n1, psi);
  //   REQUIRE(isapprox(n1c, 1.0));

  //   auto n2 = Op("Ntot", 2);
  //   complex n2c = inner(n2, psi);
  //   REQUIRE(isapprox(n2c, 2.));

  //   auto n3 = Op("Ntot", 3);
  //   complex n3c = inner(n3, psi);
  //   REQUIRE(isapprox(n3c, 0.));
  // }
  // {
  //   auto pstate = ProductState({"Up", "Dn", "Up", "Emp"});
  //   auto block = tJ(4, 2, 1);
  //   auto psi = StateReal(block, pstate);

  //   auto sz0 = Op("Sz", 0);
  //   complex sz0c = inner(sz0, psi);
  //   REQUIRE(isapprox(sz0c, 0.5));

  //   auto sz1 = Op("Sz", 1);
  //   complex sz1c = inner(sz1, psi);
  //   REQUIRE(isapprox(sz1c, -0.5));

  //   auto sz2 = Op("Sz", 2);
  //   complex sz2c = inner(sz2, psi);
  //   REQUIRE(isapprox(sz2c, 0.5));

  //   auto sz3 = Op("Sz", 3);
  //   complex sz3c = inner(sz3, psi);
  //   REQUIRE(isapprox(sz3c, 0.));

  //   auto n0 = Op("Ntot", 0);
  //   complex n0c = inner(n0, psi);
  //   REQUIRE(isapprox(n0c, 1.0));

  //   auto n1 = Op("Ntot", 1);
  //   complex n1c = inner(n1, psi);
  //   REQUIRE(isapprox(n1c, 1.0));

  //   auto n2 = Op("Ntot", 2);
  //   complex n2c = inner(n2, psi);
  //   REQUIRE(isapprox(n2c, 1.));

  //   auto n3 = Op("Ntot", 3);
  //   complex n3c = inner(n3, psi);
  //   REQUIRE(isapprox(n3c, 0.));
  // }

  // {
  //   auto pstate = ProductState({"Up", "Dn", "Up", "Dn"});
  //   auto block = Spinhalf(4, 2);
  //   auto psi = StateReal(block, pstate);

  //   auto sz0 = Op("Sz", 0);
  //   complex sz0c = inner(sz0, psi);
  //   REQUIRE(isapprox(sz0c, 0.5));

  //   auto sz1 = Op("Sz", 1);
  //   complex sz1c = inner(sz1, psi);
  //   REQUIRE(isapprox(sz1c, -0.5));

  //   auto sz2 = Op("Sz", 2);
  //   complex sz2c = inner(sz2, psi);
  //   REQUIRE(isapprox(sz2c, 0.5));

  //   auto sz3 = Op("Sz", 3);
  //   complex sz3c = inner(sz3, psi);
  //   REQUIRE(isapprox(sz3c, -0.5));
  // }
} catch (xdiag::Error const &e) {
  error_trace(e);
}
