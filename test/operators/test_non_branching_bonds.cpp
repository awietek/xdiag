#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

TEST_CASE("non_branching_bonds", "[operators]") {
  using namespace hydra;

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

  arma::cx_mat ones(arma::mat({{1.0, 1.0}, {1.0, 1.0}}),
                    arma::mat({{1.0, 1.0}, {1.0, 1.0}}));

  for (auto ss : {ones}) {
    auto bond = Bond(ss, 0);
    auto block = Spinhalf(1);
    auto h = matrix(bond, block);
    REQUIRE(close(h, ss));
  }

  // Check ground state energy of TFI model
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
  double e = eig0(bonds, block);
  double e_dmrg = -7.411918598647893;
  REQUIRE(std::abs(e - e_dmrg) < 1e-8);

  // Log.set_verbosity(2);
  // auto [ee, gs] = groundstate(bonds, block);
  // // double tau = -0.1;
  // complex tau = complex(0, -0.1);

  // auto gs2 = exp_sym_v(bonds, gs, tau);
  // HydraPrint(norm(gs2));

  // gs2.vector() /= norm(gs2);
  // complex eee = inner(bonds, gs2);
  // HydraPrint(e);
  // HydraPrint(ee);
  // HydraPrint(eee);
  // HydraPrint(norm(gs2));

  // auto allup = zero_state(block);
  // allup.vector()(block.size() - 1) = 1.0;
  // HydraPrint(norm(allup));

  // BondList mag;
  // for (int i = 0; i < N; ++i) {
  //   mag << Bond("SZ", i);
  // }

  // auto m = inner(mag, allup);
  // HydraPrint(m);

  // Log.set_verbosity(0);
  // complex delta_tau = complex(0, -0.1);
  // auto v = allup;
  // for (int t=0; t<100; ++t){
  //   v = exp_sym_v(bonds, v, delta_tau);
  //   auto m = inner(mag, v);
  //   HydraPrint(m);
  // }

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
