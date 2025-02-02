#include "../catch.hpp"

#include <iostream>
#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/operators/logic/non_branching_op.hpp>
#include <xdiag/algebra/isapprox.hpp>

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

  for (auto ss : {sx, sy, sz, sp, sm, ones}) {
    auto op = Op("Matrix", 0, ss);
    auto block = Spinhalf(1);
    auto h = matrixC(op, block);
    REQUIRE(isapprox(h, ss));
  }

  // Check ground state energy of TFI model
  int64_t N = 14;
  double J = 1.0;
  double H = 1.0;
  OpSum ops;

  for (int64_t i = 0; i < N - 1; ++i) {
    ops += J * Op("SzSz", {i, (i + 1) % N});
  }
  for (int64_t i = 0; i < N; ++i) {
    ops += H * Op("Matrix", i, sx);
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

      auto opr = Op("Matrix", sites, mr);
      auto hr = matrix(opr, block);
      REQUIRE(isapprox(hr, mr));

      auto mc = cx_mat(p2, p2, fill::randn);
      auto opc = Op("Matrix", sites, mc);
      auto hc = matrixC(opc, block);
      REQUIRE(isapprox(hc, mc));
    }
  }

  // compare scalar chirality
  {
    auto block6 = Spinhalf(6);

    OpSum ops1;
    ops1 += "Jchi" * Op("ScalarChirality", {0, 1, 2});
    ops1 += "Jchi" * Op("ScalarChirality", {1, 2, 3});
    ops1 += "Jchi" * Op("ScalarChirality", {2, 3, 4});
    ops1 += "Jchi" * Op("ScalarChirality", {3, 4, 5});
    ops1["Jchi"] = 1.0;
    auto H1 = matrixC(ops1, block6);

    cx_mat jchimat = kron(sx, kron(sy, sz) - kron(sz, sy)) +
                     kron(sy, kron(sz, sx) - kron(sx, sz)) +
                     kron(sz, kron(sx, sy) - kron(sy, sx));

    OpSum ops2;
    ops2 += Op("Matrix", {0, 1, 2}, jchimat);
    ops2 += Op("Matrix", {1, 2, 3}, jchimat);
    ops2 += Op("Matrix", {2, 3, 4}, jchimat);
    ops2 += Op("Matrix", {3, 4, 5}, jchimat);
    auto H2 = matrixC(ops2, block6);

    REQUIRE(isapprox(H1, H2));
  }

  // compare scalar chirality (12 site kagome)
  auto block12 = Spinhalf(12);

  OpSum ops1;
  ops1 += "Jchi" * Op("ScalarChirality", {0, 4, 6});
  ops1 += "Jchi" * Op("ScalarChirality", {3, 1, 9});
  ops1 += "Jchi" * Op("ScalarChirality", {9, 7, 4});
  ops1 += "Jchi" * Op("ScalarChirality", {4, 2, 10});
  ops1 += "Jchi" * Op("ScalarChirality", {10, 8, 5});
  ops1 += "Jchi" * Op("ScalarChirality", {6, 10, 1});
  ops1 += "Jchi" * Op("ScalarChirality", {1, 5, 7});
  ops1 += "Jchi" * Op("ScalarChirality", {7, 11, 2});
  ops1 += "Jchi" * Op("ScalarChirality", {2, 3, 8});
  ops1 += "Jchi" * Op("ScalarChirality", {8, 9, 0});
  ops1 += "Jchi" * Op("ScalarChirality", {5, 0, 11});
  ops1 += "Jchi" * Op("ScalarChirality", {11, 6, 3});
  ops1 += "Jchi" * Op("ScalarChirality", {4, 10, 6});
  ops1 += "Jchi" * Op("ScalarChirality", {1, 7, 9});
  ops1 += "Jchi" * Op("ScalarChirality", {7, 2, 4});
  ops1 += "Jchi" * Op("ScalarChirality", {2, 8, 10});
  ops1 += "Jchi" * Op("ScalarChirality", {8, 0, 5});
  ops1 += "Jchi" * Op("ScalarChirality", {10, 5, 1});
  ops1 += "Jchi" * Op("ScalarChirality", {5, 11, 7});
  ops1 += "Jchi" * Op("ScalarChirality", {11, 3, 2});
  ops1 += "Jchi" * Op("ScalarChirality", {3, 9, 8});
  ops1 += "Jchi" * Op("ScalarChirality", {9, 4, 0});
  ops1 += "Jchi" * Op("ScalarChirality", {0, 6, 11});
  ops1 += "Jchi" * Op("ScalarChirality", {6, 1, 3});
  ops1["Jchi"] = 1.0;
  auto H1 = matrixC(ops1, block12);

  cx_mat jchimat = kron(kron(sx, sy), sz) + kron(kron(sy, sz), sx) +
                   kron(kron(sz, sx), sy) - kron(kron(sx, sz), sy) -
                   kron(kron(sy, sx), sz) - kron(kron(sz, sy), sx);

  OpSum ops2;
  ops2 += Op("Matrix", {0, 4, 6}, jchimat);
  ops2 += Op("Matrix", {3, 1, 9}, jchimat);
  ops2 += Op("Matrix", {9, 7, 4}, jchimat);
  ops2 += Op("Matrix", {4, 2, 10}, jchimat);
  ops2 += Op("Matrix", {10, 8, 5}, jchimat);
  ops2 += Op("Matrix", {6, 10, 1}, jchimat);
  ops2 += Op("Matrix", {1, 5, 7}, jchimat);
  ops2 += Op("Matrix", {7, 11, 2}, jchimat);
  ops2 += Op("Matrix", {2, 3, 8}, jchimat);
  ops2 += Op("Matrix", {8, 9, 0}, jchimat);
  ops2 += Op("Matrix", {5, 0, 11}, jchimat);
  ops2 += Op("Matrix", {11, 6, 3}, jchimat);
  ops2 += Op("Matrix", {4, 10, 6}, jchimat);
  ops2 += Op("Matrix", {1, 7, 9}, jchimat);
  ops2 += Op("Matrix", {7, 2, 4}, jchimat);
  ops2 += Op("Matrix", {2, 8, 10}, jchimat);
  ops2 += Op("Matrix", {8, 0, 5}, jchimat);
  ops2 += Op("Matrix", {10, 5, 1}, jchimat);
  ops2 += Op("Matrix", {5, 11, 7}, jchimat);
  ops2 += Op("Matrix", {11, 3, 2}, jchimat);
  ops2 += Op("Matrix", {3, 9, 8}, jchimat);
  ops2 += Op("Matrix", {9, 4, 0}, jchimat);
  ops2 += Op("Matrix", {0, 6, 11}, jchimat);
  ops2 += Op("Matrix", {6, 1, 3}, jchimat);
  auto H2 = matrixC(ops2, block12);

  // XDIAG_SHOW(norm(H1));
  // XDIAG_SHOW(norm(H2));

  // REQUIRE(isapprox(H1, H2));
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
