#include "../catch.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/utils/close.hpp>

using namespace xdiag;

TEST_CASE("product_state_distributed", "[states]") {

  // Test product state for
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
      // XDIAG_SHOW(pstate);

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

      auto block2 = tJDistributed(n_sites, nup, ndn);
      auto psi2 = State(block2);
      fill(psi2, pstate);
      for (int i = 0; i < n_sites; ++i) {
        std::string p = pstate[i];
        auto sz = Op("SZ", i);
        auto szc = inner(sz, psi2);

        auto n = Op("NTOT", i);
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
      // XDIAG_SHOW(pstate);

      int nup = 0;
      auto block = Spinhalf(n_sites);
      auto psi = State(block);
      fill(psi, pstate);
      for (int i = 0; i < n_sites; ++i) {
        std::string p = pstate[i];

        auto sz = Op("SZ", i);
        auto szc = inner(sz, psi);

        if (p == "Up") {
          REQUIRE(close(szc, 0.5));
          ++nup;
        } else if (p == "Dn") {
          REQUIRE(close(szc, -0.5));
        }
      }

      auto block2 = SpinhalfDistributed(n_sites, nup);
      auto psi2 = State(block2);
      fill(psi2, pstate);
      for (int i = 0; i < n_sites; ++i) {
        std::string p = pstate[i];
        auto sz = Op("SZ", i);
        auto szc = inner(sz, psi2);

        if (p == "Up") {
          REQUIRE(close(szc, 0.5));
        } else if (p == "Dn") {
          REQUIRE(close(szc, -0.5));
        }
      }
    }
  }
}
