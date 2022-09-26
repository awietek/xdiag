#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

TEST_CASE("time_evolve", "[algorithms]") {
  using namespace hydra;
  using namespace arma;

  Log("time_evolve test");
  int N = 10;
  double J = 1.0;
  double H = 1.0;
  BondList bonds;

  cx_mat sx(mat({{0., 0.5}, {0.5, 0.}}), mat({{0., 0.}, {0., 0.}}));

  for (int i = 0; i < N - 1; ++i) {
    bonds << Bond("ISING", J, {i, (i + 1) % N});
  }
  for (int i = 0; i < N; ++i) {
    bonds << Bond(sx, H, i);
  }

  auto block = Spinhalf(N);
  auto Hmat = matrix_real(bonds, block);

  double tau = 0.1;

  auto v0 = zero_state(block);
  v0.vector()(block.size() - 1) = 1.0;
  auto v = time_evolve(bonds, v0, tau);
  arma::cx_mat eHv = expmat(complex(0, -tau) * Hmat) * v0.vector();

  HydraPrint(eHv);
  HydraPrint(v.vector());
  HydraPrint(norm(eHv - v.vector()));
  REQUIRE(norm(eHv - v.vector()) < 1e-8);
}
