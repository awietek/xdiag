// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>

TEST_CASE("electron_raise_lower", "[electron]") try {
  using namespace xdiag;

  // Test normal ordering
  for (int nsites = 2; nsites < 6; ++nsites) {
    auto block0 = Electron(nsites, 0, 0);
    auto psi0 = State(block0);

    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block = Electron(nsites, nup, ndn);

        for (auto pstate : block) {
          std::vector<int> up_positions;
          std::vector<int> dn_positions;
          for (int i = 0; i < nsites; ++i) {
            if ((pstate[i] == "Up") || (pstate[i] == "UpDn")) {
              up_postitions.push_back(i);
            } else if ((pstate[i] == "Dn") || (pstate[i] == "UpDn")) {
              dn_postitions.push_back(i);
            }
          }

          // Create state from product state
          auto psi = State(block);
          fill(psi, pstate);

          auto psi2 = psi0;
          for (int i :
               std::vector<int>(up_positions.rbegin(), up_positions.rend())) {
            psi2 = apply(Op("Cdagup", i), psi2);
          }
          for (int i :
               std::vector<int>(dn_positions.rbegin(), dn_positions.rend())) {
            psi2 = apply(Op("Cdagdn", i), psi2);
          }
	  REQUIRE(isclose(psi2, psi));
	  
        }
      }
    }
  }

  using namespace arma;

  std::vector<std::string> op_strs = {"Cdagup", "Cdagdn", "Cup", "Cdn"};

  for (int nsites = 1; nsites < 4; ++nsites) {
    Log("testing electron anticommutation relations (full matrix): N={}",
        nsites);

    auto block = Electron(nsites);
    int64_t D = block.size();

    cx_mat id(D, D, fill::eye);

    for (int i = 0; i < nsites; ++i) {
      for (int j = 0; j < nsites; ++j) {

        for (auto op_i_str : op_strs) {
          for (auto op_j_str : op_strs) {
            auto op_i = Op(op_i_str, i);
            auto op_j = Op(op_j_str, j);
            auto op_i_m = matrixC(op_i, block);
            auto op_j_m = matrixC(op_j, block);
            cx_mat anti_comm = op_i_m * op_j_m + op_j_m * op_i_m;

            // Log("{} {}; {} {}", op_i_str, i, op_j_str, j);

            if (((op_i_str == "Cdagup") && (op_j_str == "Cup")) ||
                ((op_i_str == "Cdagdn") && (op_j_str == "Cdn")) ||
                ((op_i_str == "Cup") && (op_j_str == "Cdagup")) ||
                ((op_i_str == "Cdn") && (op_j_str == "Cdagdn"))) {
              if (i == j) {
                REQUIRE(norm(anti_comm - id) < 1e-12);
              } else {
                REQUIRE(norm(anti_comm) < 1e-12);
              }
            } else {
              REQUIRE(norm(anti_comm) < 1e-12);
            }
          }
        }
      }
    }
  }

  for (int nsites = 2; nsites < 4; ++nsites) {

    auto block = Electron(nsites);
    int64_t D = block.size();

    for (int i = 0; i < nsites; ++i) {

      // Check whether Nup agrees
      {
        auto op = Op("Nup", i);
        auto m = matrix(op, block);
        auto mi1 = matrix(Op("Cdagup", i), block);
        auto mi2 = matrix(Op("Cup", i), block);
        auto mm = mi1 * mi2;
        REQUIRE(norm(mm - m) < 1e-12);
      }

      // Check whether Ndn agrees
      {
        auto op = Op("Ndn", i);
        auto m = matrix(op, block);
        auto mi1 = matrix(Op("Cdagdn", i), block);
        auto mi2 = matrix(Op("Cdn", i), block);
        auto mm = mi1 * mi2;
        REQUIRE(norm(mm - m) < 1e-12);
      }

      for (int j = 0; j < nsites; ++j) {
        if (i != j) {

          // Check whether Hopup agrees
          {
            auto op = Op("Hopup", {i, j});
            auto m = matrix(op, block);
            auto m1 = matrix(Op("Cdagup", i), block);
            auto m2 = matrix(Op("Cup", j), block);
            auto m3 = matrix(Op("Cdagup", j), block);
            auto m4 = matrix(Op("Cup", i), block);
            auto mm = -m1 * m2 - m3 * m4;
            REQUIRE(norm(mm - m) < 1e-12);
          }

          // Check whether Hopdn agrees
          {
            auto op = Op("Hopdn", {i, j});
            auto m = matrix(op, block);
            auto m1 = matrix(Op("Cdagdn", i), block);
            auto m2 = matrix(Op("Cdn", j), block);
            auto m3 = matrix(Op("Cdagdn", j), block);
            auto m4 = matrix(Op("Cdn", i), block);
            auto mm = -m1 * m2 - m3 * m4;
            REQUIRE(norm(mm - m) < 1e-12);
          }

          // Check whether SzSz agrees
          {
            auto op = Op("SzSz", {i, j});
            auto m = matrix(op, block);

            auto mi1 = matrix(Op("Cdagup", i), block);
            auto mi2 = matrix(Op("Cup", i), block);
            auto mi3 = matrix(Op("Cdagdn", i), block);
            auto mi4 = matrix(Op("Cdn", i), block);
            auto szi = 0.5 * mi1 * mi2 - 0.5 * mi3 * mi4;

            auto mj1 = matrix(Op("Cdagup", j), block);
            auto mj2 = matrix(Op("Cup", j), block);
            auto mj3 = matrix(Op("Cdagdn", j), block);
            auto mj4 = matrix(Op("Cdn", j), block);
            auto szj = 0.5 * mj1 * mj2 - 0.5 * mj3 * mj4;

            auto mm = szi * szj;
            REQUIRE(norm(mm - m) < 1e-12);
          }

          // Check whether Exchange agrees
          {
            auto op = Op("Exchange", {i, j});
            auto m = matrix(op, block);

            auto mi1 = matrix(Op("Cdagup", i), block);
            auto mi2 = matrix(Op("Cup", i), block);
            auto mi3 = matrix(Op("Cdagdn", i), block);
            auto mi4 = matrix(Op("Cdn", i), block);
            auto spi = mi1 * mi4;
            auto smi = mi3 * mi2;

            auto mj1 = matrix(Op("Cdagup", j), block);
            auto mj2 = matrix(Op("Cup", j), block);
            auto mj3 = matrix(Op("Cdagdn", j), block);
            auto mj4 = matrix(Op("Cdn", j), block);
            auto spj = mj1 * mj4;
            auto smj = mj3 * mj2;

            auto mm = 0.5 * spi * smj + 0.5 * smi * spj;
            REQUIRE(norm(mm - m) < 1e-12);
          }
        }
      }
    }
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
}
