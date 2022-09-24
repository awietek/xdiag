#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

TEST_CASE("non_branching_bonds", "[operators]") {
  using namespace hydra;

  // Log("hello");
  arma::cx_mat sx(arma::mat({{0., 0.5}, {0.5, 0.}}),
                  arma::mat({{0., 0.}, {0., 0.}}));
  arma::cx_mat sy(arma::mat({{0., 0.}, {0., 0.}}),
                  arma::mat({{0., -0.5}, {0.5, 0.}}));
  arma::cx_mat sz(arma::mat({{0.5, 0.0}, {0.0, 0.5}}),
                  arma::mat({{0., 0.}, {0., 0.0}}));

  arma::cx_mat sp(arma::mat({{0.0, 1.0}, {0.0, 0.0}}),
                  arma::mat({{0., 0.}, {0., 0.0}}));

  arma::cx_mat sm(arma::mat({{0.0, 0.0}, {1.0, 0.0}}),
                  arma::mat({{0., 0.}, {0., 0.0}}));

  for (auto ss : {sx, sy, sz, sp, sm}) {
    auto bond = Bond(ss, 0);
    auto block = Spinhalf(1);
    auto h = matrix(bond, block);
    REQUIRE(close(h, ss));
  }

  int N = 14;
  double J = 1.0;
  double H = 1.0;
  BondList bonds;

  for (int i = 0; i < N - 1; ++i) {
    bonds << Bond("ISING", J, {i, (i + 1) % N});
  }
  for (int i = 0; i < N; ++i) {
    bonds << Bond(sx, H, i);
  }

  auto block = Spinhalf(N);
  double e = e0(bonds, block);
  double e_dmrg = -7.411918598647893;
  // HydraPrint(e);
  // HydraPrint(e_dmrg);
  REQUIRE(std::abs(e - e_dmrg) < 1e-8);

  // HydraPrint(sx);
  // HydraPrint(sy);
  // HydraPrint(sz);

  // auto sx_bond = Bond(sx, 0);
  // REQUIRE(operators::is_non_branching_bond(sx_bond));

  // auto block = Spinhalf(1);
  // HydraPrint(result);

  // arma::mat sxsx = arma::kron(sx, sx);
  // auto sxsx_bond = Bond(sxsx, {0, 1});
  // auto block2 = Spinhalf(2);
  // auto result2 = matrix_real(sxsx_bond, block2);
  // HydraPrint(sxsx);
  // HydraPrint(result2);

  // arma::mat sxsx = arma::kron(sx, sx);
  // auto sxsx_bond = Bond(sxsx, {0, 1});
  // auto block2 = Spinhalf(2);
  // auto result2 = matrix_real(sxsx_bond, block2);
  // HydraPrint(sxsx);
  // HydraPrint(result2);
}
