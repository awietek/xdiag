#include "../catch.hpp"

#include <iostream>

#include "testcases_tj.h"
#include <hydra/all.h>

using namespace hydra;

void test_tjmodel_e0_real(BondList bonds, Couplings couplings, int nup, int ndn,
                          double e0) {
  int n_sites = bonds.n_sites();

  auto block = tJ<uint32>(n_sites, nup, ndn);
  auto H = matrix_real(bonds, couplings, block, block);
  REQUIRE(lila::close(H, lila::Herm(H)));

  auto eigs = lila::EigenvaluesSym(H);
  HydraLog.out("{} {} {} {}", nup, ndn, eigs(0), e0);
  CHECK(std::abs(e0 - eigs(0)) < 1e-6);
}

TEST_CASE("tj_matrix", "[models]") {
  using namespace hydra::testcases::tj;

  {
    HydraLog.out("TJModel: six-site chain test, t=1.0, J=1.0, N=6");
    auto [bonds, cpls] = tJchain(6, 1.0, 1.0);
    std::vector<std::tuple<int, int, double>> nup_ndn_e0 = {
        {0, 0, 0.0},         {0, 1, -2.0},        {0, 2, -2.96081311},
        {0, 3, -3.79610527}, {0, 4, -2.46081311}, {0, 5, -0.99999999},
        {0, 6, 1.500000000}, {1, 0, -2.0},        {1, 1, -3.61222054},
        {1, 2, -4.04537829}, {1, 3, -4.10768318}, {1, 4, -2.42705097},
        {1, 5, -0.49999999}, {2, 0, -2.96081311}, {2, 1, -4.04537829},
        {2, 2, -4.16447847}, {2, 3, -3.52922048}, {2, 4, -2.11803398},
        {3, 0, -3.79610527}, {3, 1, -4.10768318}, {3, 2, -3.52922048},
        {3, 3, -2.80277563}, {4, 0, -2.46081311}, {4, 1, -2.42705097},
        {4, 2, -2.11803398}, {5, 0, -0.99999999}, {5, 1, -0.49999999},
        {6, 0, 1.500000000}};
    for (auto [nup, ndn, e0] : nup_ndn_e0)
      test_tjmodel_e0_real(bonds, cpls, nup, ndn, e0);
  }

  {
    HydraLog.out("TJModel: six-site chain test, t=1.0, J=0.0, N=6");
    auto [bonds, cpls] = tJchain(6, 1.0, 0.0);
    std::vector<std::tuple<int, int, double>> nup_ndn_e0 = {
        {0, 0, 0.0},         {0, 1, -2.0},        {0, 2, -3.00000000},
        {0, 3, -4.00000000}, {0, 4, -2.99999999}, {0, 5, -2.00000000},
        {0, 6, 0.000000000}, {1, 0, -2.0},        {1, 1, -3.46410161},
        {1, 2, -3.99999999}, {1, 3, -3.46410161}, {1, 4, -1.99999999},
        {1, 5, 0.000000000}, {2, 0, -3.00000000}, {2, 1, -3.99999999},
        {2, 2, -3.46410161}, {2, 3, -1.99999999}, {2, 4, 0.000000000},
        {3, 0, -4.00000000}, {3, 1, -3.46410161}, {3, 2, -1.99999999},
        {3, 3, 0.000000000}, {4, 0, -2.99999999}, {4, 1, -1.99999999},
        {4, 2, 0.000000000}, {5, 0, -2.00000000}, {5, 1, 0.000000000},
        {6, 0, 0.000000000}};
    for (auto [nup, ndn, e0] : nup_ndn_e0)
      test_tjmodel_e0_real(bonds, cpls, nup, ndn, e0);
  }
}
