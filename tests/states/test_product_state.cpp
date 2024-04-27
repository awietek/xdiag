#include "../catch.hpp"

#include <xdiag/states/product_state.hpp>
#include <xdiag/algebra/algebra.hpp>
#include <xdiag/utils/close.hpp>

using namespace xdiag;

TEST_CASE("product_state", "[states]") {

  // Test product state for Electron block
  for (int i = 0; i < 10; ++i) {
    for (int n_sites = 1; n_sites < 9; ++n_sites) {

      std::vector<std::string> ps = {"Emp", "Up", "Dn", "UpDn",
                                     "Emp", "Up", "Dn", "UpDn"};
      std::vector<std::string> pss;
      std::sample(ps.begin(), ps.end(), std::back_inserter(pss), n_sites,
                  std::mt19937(i));
      auto pstate = ProductState(pss);
      // XDIAG_PRINT(pstate);

      auto block = Electron(n_sites);
      auto psi = State(block, true);
      fill(psi, pstate);
      
      int nup = 0;
      int ndn = 0;
      for (int i = 0; i < n_sites; ++i) {
        std::string p = pstate[i];
        auto sz = Bond("SZ", i);
        auto szc = inner(sz, psi);

        auto n = Bond("NUMBER", i);
        auto nc = inner(n, psi);

        if (p == "Emp") {
          REQUIRE(close(szc, 0.));
          REQUIRE(close(nc, 0.));
        } else if (p == "Up") {
          REQUIRE(close(szc, 0.5));
          REQUIRE(close(nc, 1.0));
          ++nup;
        } else if (p == "Dn") {
          REQUIRE(close(szc, -0.5));
          REQUIRE(close(nc, 1.0));
          ++ndn;
        } else if (p == "UpDn") {
          REQUIRE(close(szc, 0.0));
          REQUIRE(close(nc, 2.0));
          ++nup;
          ++ndn;
        }
      }

      auto block2 = Electron(n_sites, nup, ndn);
      auto psi2 = State(block2);
      fill(psi2, pstate);
      
      for (int i = 0; i < n_sites; ++i) {
        std::string p = pstate[i];
        auto sz = Bond("SZ", i);
        auto szc = inner(sz, psi2);

        auto n = Bond("NUMBER", i);
        auto nc = inner(n, psi2);

        if (p == "Emp") {
          REQUIRE(close(szc, 0.));
          REQUIRE(close(nc, 0.));
        } else if (p == "Up") {
          REQUIRE(close(szc, 0.5));
          REQUIRE(close(nc, 1.0));
        } else if (p == "Dn") {
          REQUIRE(close(szc, -0.5));
          REQUIRE(close(nc, 1.0));
        } else if (p == "UpDn") {
          REQUIRE(close(szc, 0.0));
          REQUIRE(close(nc, 2.0));
        }
      }
    }
  }

  // Test product state for tJ block
  for (int i = 0; i < 10; ++i) {
    for (int n_sites = 1; n_sites < 9; ++n_sites) {

      std::vector<std::string> ps = {
          "Emp", "Up", "Dn", "Emp", "Up", "Dn", "Emp", "Up", "Dn",
      };
      std::vector<std::string> pss;
      std::sample(ps.begin(), ps.end(), std::back_inserter(pss), n_sites,
                  std::mt19937(i));
      auto pstate = ProductState(pss);
      // XDIAG_PRINT(pstate);

      int nup = 0;
      int ndn = 0;
      for (int i = 0; i < n_sites; ++i) {
        std::string p = pstate[i];

        if (p == "Up") {
          ++nup;
        } else if (p == "Dn") {
          ++ndn;
        }
      }

      auto block2 = tJ(n_sites, nup, ndn);
      auto psi2 = State(block2);
      fill(psi2, pstate);
      for (int i = 0; i < n_sites; ++i) {
        std::string p = pstate[i];
        auto sz = Bond("SZ", i);
        auto szc = inner(sz, psi2);

        auto n = Bond("NUMBER", i);
        auto nc = inner(n, psi2);

        if (p == "Emp") {
          REQUIRE(close(szc, 0.));
          REQUIRE(close(nc, 0.));
        } else if (p == "Up") {
          REQUIRE(close(szc, 0.5));
          REQUIRE(close(nc, 1.0));
        } else if (p == "Dn") {
          REQUIRE(close(szc, -0.5));
          REQUIRE(close(nc, 1.0));
        }
      }
    }
  }

  // Test product state for Spinhalf block
  for (int i = 0; i < 10; ++i) {
    for (int n_sites = 1; n_sites < 9; ++n_sites) {

      std::vector<std::string> ps = {"Up", "Dn", "Up", "Dn",
                                     "Up", "Dn", "Up", "Dn"};
      std::vector<std::string> pss;
      std::sample(ps.begin(), ps.end(), std::back_inserter(pss), n_sites,
                  std::mt19937(i));
      auto pstate = ProductState(pss);
      // XDIAG_PRINT(pstate);

      int nup = 0;
      auto block = Spinhalf(n_sites);
      auto psi = State(block);
      fill(psi, pstate);
      for (int i = 0; i < n_sites; ++i) {
        std::string p = pstate[i];

        auto sz = Bond("SZ", i);
        auto szc = inner(sz, psi);
	
        if (p == "Up") {
          REQUIRE(close(szc, 0.5));
          ++nup;
        } else if (p == "Dn") {
	  REQUIRE(close(szc, -0.5));
        }
      }

      auto block2 = Spinhalf(n_sites, nup);
      auto psi2 = State(block2);
      fill(psi2, pstate);
      for (int i = 0; i < n_sites; ++i) {
        std::string p = pstate[i];
        auto sz = Bond("SZ", i);
        auto szc = inner(sz, psi2);

        if (p == "Up") {
          REQUIRE(close(szc, 0.5));
        } else if (p == "Dn") {
          REQUIRE(close(szc, -0.5));
        }
      }
    }
  }

  // {
  //   for

  //     auto pstate = ProductState({"Up", "Dn", "UpDn", "Emp"});
  //   auto block = Electron(4, 2, 2);
  //   auto psi = StateReal(block, pstate);

  //   auto sz0 = Bond("SZ", 0);
  //   complex sz0c = inner(sz0, psi);
  //   REQUIRE(close(sz0c, 0.5));

  //   auto sz1 = Bond("SZ", 1);
  //   complex sz1c = inner(sz1, psi);
  //   REQUIRE(close(sz1c, -0.5));

  //   auto sz2 = Bond("SZ", 2);
  //   complex sz2c = inner(sz2, psi);
  //   REQUIRE(close(sz2c, 0.));

  //   auto sz3 = Bond("SZ", 3);
  //   complex sz3c = inner(sz3, psi);
  //   REQUIRE(close(sz3c, 0.));

  //   auto n0 = Bond("NUMBER", 0);
  //   complex n0c = inner(n0, psi);
  //   REQUIRE(close(n0c, 1.0));

  //   auto n1 = Bond("NUMBER", 1);
  //   complex n1c = inner(n1, psi);
  //   REQUIRE(close(n1c, 1.0));

  //   auto n2 = Bond("NUMBER", 2);
  //   complex n2c = inner(n2, psi);
  //   REQUIRE(close(n2c, 2.));

  //   auto n3 = Bond("NUMBER", 3);
  //   complex n3c = inner(n3, psi);
  //   REQUIRE(close(n3c, 0.));
  // }
  // {
  //   auto pstate = ProductState({"Up", "Dn", "Up", "Emp"});
  //   auto block = tJ(4, 2, 1);
  //   auto psi = StateReal(block, pstate);

  //   auto sz0 = Bond("SZ", 0);
  //   complex sz0c = inner(sz0, psi);
  //   REQUIRE(close(sz0c, 0.5));

  //   auto sz1 = Bond("SZ", 1);
  //   complex sz1c = inner(sz1, psi);
  //   REQUIRE(close(sz1c, -0.5));

  //   auto sz2 = Bond("SZ", 2);
  //   complex sz2c = inner(sz2, psi);
  //   REQUIRE(close(sz2c, 0.5));

  //   auto sz3 = Bond("SZ", 3);
  //   complex sz3c = inner(sz3, psi);
  //   REQUIRE(close(sz3c, 0.));

  //   auto n0 = Bond("NUMBER", 0);
  //   complex n0c = inner(n0, psi);
  //   REQUIRE(close(n0c, 1.0));

  //   auto n1 = Bond("NUMBER", 1);
  //   complex n1c = inner(n1, psi);
  //   REQUIRE(close(n1c, 1.0));

  //   auto n2 = Bond("NUMBER", 2);
  //   complex n2c = inner(n2, psi);
  //   REQUIRE(close(n2c, 1.));

  //   auto n3 = Bond("NUMBER", 3);
  //   complex n3c = inner(n3, psi);
  //   REQUIRE(close(n3c, 0.));
  // }

  // {
  //   auto pstate = ProductState({"Up", "Dn", "Up", "Dn"});
  //   auto block = Spinhalf(4, 2);
  //   auto psi = StateReal(block, pstate);

  //   auto sz0 = Bond("SZ", 0);
  //   complex sz0c = inner(sz0, psi);
  //   REQUIRE(close(sz0c, 0.5));

  //   auto sz1 = Bond("SZ", 1);
  //   complex sz1c = inner(sz1, psi);
  //   REQUIRE(close(sz1c, -0.5));

  //   auto sz2 = Bond("SZ", 2);
  //   complex sz2c = inner(sz2, psi);
  //   REQUIRE(close(sz2c, 0.5));

  //   auto sz3 = Bond("SZ", 3);
  //   complex sz3c = inner(sz3, psi);
  //   REQUIRE(close(sz3c, -0.5));
  // }
}
