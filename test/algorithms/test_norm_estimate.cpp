#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

#include "../blocks/spinhalf/testcases_spinhalf.h"


TEST_CASE("norm_estimate", "[algorithms]") {
  using namespace hydra;
  using namespace arma;
  int n = 200;
  Log.set_verbosity(2);

  for (int seed = 1; seed <= 10; ++seed) {
    arma_rng::set_seed(seed);
    mat B(n, n, fill::randn);
    mat A = B + B.t();
    double norm_exact = norm(A, "inf");
    double norm_est =
        norm_estimate([&A](arma::vec const &v) { return mat(A * v); }, n);
    double ratio = norm_exact / norm_est;
    REQUIRE(((ratio < 2) && (ratio > 0.5)));

    Log("norm_exact: {}", norm_exact);
    Log("norm_est: {}", norm_est);
  }
  int n_sites = 10;
  auto bonds = hydra::testcases::spinhalf::HB_alltoall(n_sites);
  auto block = Spinhalf(n_sites, 5);
  auto H = matrix_real(bonds, block, block);
  double norm_exact = norm(H, "inf");
  double norm_est =
    norm_estimate([&H](arma::vec const &v) { return mat(H * v); }, block.size());
  double ratio = norm_exact / norm_est;
  REQUIRE(((ratio < 2) && (ratio > 0.5)));

  Log("norm_exact: {}", norm_exact);
  Log("norm_est: {}", norm_est);
}
