#include "../catch.hpp"

#include <iostream>
#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/operators/non_branching_op.hpp>
#include <xdiag/utils/close.hpp>

TEST_CASE("non_branching_op", "[operators]") try {
  using namespace xdiag;
  using namespace arma;
  Log("testing non_branching_op");

  cx_mat sx(mat({{0., 0.5}, {0.5, 0.}}), mat({{0., 0.}, {0., 0.}}));
  cx_mat sy(mat({{0., 0.}, {0., 0.}}), mat({{0., -0.5}, {0.5, 0.}}));
  cx_mat sz(mat({{0.5, 0.0}, {0.0, -0.5}}), mat({{0., 0.}, {0., 0.0}}));

  cx_mat sp(mat({{0.0, 1.0}, {0.0, 0.0}}), mat({{0., 0.}, {0., 0.0}}));
  cx_mat sm(mat({{0.0, 0.0}, {1.0, 0.0}}), mat({{0., 0.}, {0., 0.0}}));
  cx_mat ones(mat({{1.0, 1.0}, {1.0, 1.0}}), mat({{1.0, 1.0}, {1.0, 1.0}}));

  for (auto ss : {ones}) {
    auto op = Op("Map", ss, 0);
    auto block = Spinhalf(1);
    auto h = matrixC(op, block);
    REQUIRE(close(h, ss));
  }

  // Check ground state energy of TFI model
  int64_t N = 14;
  double J = 1.0;
  double H = 1.0;
  OpSum ops;

  for (int64_t i = 0; i < N - 1; ++i) {
    ops += Op("SzSz", J, {i, (i + 1) % N});
  }
  for (int64_t i = 0; i < N; ++i) {
    ops += Op("SX", sx, i);
  }

  auto block = Spinhalf(N);
  double e = eigval0(ops, block);
  double e_dmrg = -7.411918598647893;
  REQUIRE(std::abs(e - e_dmrg) < 1e-8);

  // Check whether random ops are created correctly
  for (int64_t r = 0; r < 5; ++r) {

    for (int64_t k = 1; k <= 6; ++k) {
      auto block = Spinhalf(k);
      int64_t p2 = pow(2, k);

      auto mr = mat(p2, p2, fill::randn);
      std::vector<int64_t> sites(k);
      std::iota(sites.begin(), sites.end(), 0);

      auto opr = Op("MR", mr, sites);
      auto hr = matrix(opr, block);
      REQUIRE(close(hr, mr));

      auto mc = cx_mat(p2, p2, fill::randn);
      auto opc = Op("MC", mc, sites);
      auto hc = matrixC(opc, block);
      REQUIRE(close(hc, mc));
    }
  }

  // compare scalar chirality
  {
    auto block6 = Spinhalf(6);

    OpSum ops1;
    ops1 += Op("ScalarChirality", "Jchi", {0, 1, 2});
    ops1 += Op("ScalarChirality", "Jchi", {1, 2, 3});
    ops1 += Op("ScalarChirality", "Jchi", {2, 3, 4});
    ops1 += Op("ScalarChirality", "Jchi", {3, 4, 5});
    ops1["Jchi"] = 1.0;
    auto H1 = matrixC(ops1, block6);

    cx_mat jchimat = kron(sx, kron(sy, sz) - kron(sz, sy)) +
                     kron(sy, kron(sz, sx) - kron(sx, sz)) +
                     kron(sz, kron(sx, sy) - kron(sy, sx));

    OpSum ops2;
    ops2 += Op("Jchi", jchimat, {0, 1, 2});
    ops2 += Op("Jchi", jchimat, {1, 2, 3});
    ops2 += Op("Jchi", jchimat, {2, 3, 4});
    ops2 += Op("Jchi", jchimat, {3, 4, 5});
    auto H2 = matrixC(ops2, block6);

    REQUIRE(close(H1, H2));
  }

  // compare scalar chirality (12 site kagome)
  auto block12 = Spinhalf(12);

  OpSum ops1;
  ops1 += Op("ScalarChirality", "Jchi", {0, 4, 6});
  ops1 += Op("ScalarChirality", "Jchi", {3, 1, 9});
  ops1 += Op("ScalarChirality", "Jchi", {9, 7, 4});
  ops1 += Op("ScalarChirality", "Jchi", {4, 2, 10});
  ops1 += Op("ScalarChirality", "Jchi", {10, 8, 5});
  ops1 += Op("ScalarChirality", "Jchi", {6, 10, 1});
  ops1 += Op("ScalarChirality", "Jchi", {1, 5, 7});
  ops1 += Op("ScalarChirality", "Jchi", {7, 11, 2});
  ops1 += Op("ScalarChirality", "Jchi", {2, 3, 8});
  ops1 += Op("ScalarChirality", "Jchi", {8, 9, 0});
  ops1 += Op("ScalarChirality", "Jchi", {5, 0, 11});
  ops1 += Op("ScalarChirality", "Jchi", {11, 6, 3});
  ops1 += Op("ScalarChirality", "Jchi", {4, 10, 6});
  ops1 += Op("ScalarChirality", "Jchi", {1, 7, 9});
  ops1 += Op("ScalarChirality", "Jchi", {7, 2, 4});
  ops1 += Op("ScalarChirality", "Jchi", {2, 8, 10});
  ops1 += Op("ScalarChirality", "Jchi", {8, 0, 5});
  ops1 += Op("ScalarChirality", "Jchi", {10, 5, 1});
  ops1 += Op("ScalarChirality", "Jchi", {5, 11, 7});
  ops1 += Op("ScalarChirality", "Jchi", {11, 3, 2});
  ops1 += Op("ScalarChirality", "Jchi", {3, 9, 8});
  ops1 += Op("ScalarChirality", "Jchi", {9, 4, 0});
  ops1 += Op("ScalarChirality", "Jchi", {0, 6, 11});
  ops1 += Op("ScalarChirality", "Jchi", {6, 1, 3});
  ops1["Jchi"] = 1.0;
  auto H1 = matrixC(ops1, block12);

  cx_mat jchimat = kron(kron(sx, sy), sz) + kron(kron(sy, sz), sx) +
                   kron(kron(sz, sx), sy) - kron(kron(sx, sz), sy) -
                   kron(kron(sy, sx), sz) - kron(kron(sz, sy), sx);

  OpSum ops2;
  ops2 += Op("Jchi", jchimat, {0, 4, 6});
  ops2 += Op("Jchi", jchimat, {3, 1, 9});
  ops2 += Op("Jchi", jchimat, {9, 7, 4});
  ops2 += Op("Jchi", jchimat, {4, 2, 10});
  ops2 += Op("Jchi", jchimat, {10, 8, 5});
  ops2 += Op("Jchi", jchimat, {6, 10, 1});
  ops2 += Op("Jchi", jchimat, {1, 5, 7});
  ops2 += Op("Jchi", jchimat, {7, 11, 2});
  ops2 += Op("Jchi", jchimat, {2, 3, 8});
  ops2 += Op("Jchi", jchimat, {8, 9, 0});
  ops2 += Op("Jchi", jchimat, {5, 0, 11});
  ops2 += Op("Jchi", jchimat, {11, 6, 3});
  ops2 += Op("Jchi", jchimat, {4, 10, 6});
  ops2 += Op("Jchi", jchimat, {1, 7, 9});
  ops2 += Op("Jchi", jchimat, {7, 2, 4});
  ops2 += Op("Jchi", jchimat, {2, 8, 10});
  ops2 += Op("Jchi", jchimat, {8, 0, 5});
  ops2 += Op("Jchi", jchimat, {10, 5, 1});
  ops2 += Op("Jchi", jchimat, {5, 11, 7});
  ops2 += Op("Jchi", jchimat, {11, 3, 2});
  ops2 += Op("Jchi", jchimat, {3, 9, 8});
  ops2 += Op("Jchi", jchimat, {9, 4, 0});
  ops2 += Op("Jchi", jchimat, {0, 6, 11});
  ops2 += Op("Jchi", jchimat, {6, 1, 3});
  auto H2 = matrixC(ops2, block12);

  // XDIAG_SHOW(norm(H1));
  // XDIAG_SHOW(norm(H2));

  REQUIRE(close(H1, H2));
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
