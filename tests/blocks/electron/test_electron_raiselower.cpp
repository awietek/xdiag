#include "../../catch.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/blocks/electron/apply.hpp>
#include <xdiag/blocks/electron/matrix.hpp>
#include <xdiag/utils/close.hpp>

TEST_CASE("electron_raise_lower", "[electron]") {
  using namespace xdiag;
  using namespace arma;
  std::vector<std::string> op_strs = {"CDAGUP", "CDAGDN", "CUP", "CDN"};

  for (int n_sites = 1; n_sites < 4; ++n_sites) {
    Log("testing electron anticommutation relations (full matrix): N={}",
        n_sites);

    auto block = Electron(n_sites);
    int64_t D = block.size();

    cx_mat id(D, D, fill::eye);

    for (int i = 0; i < n_sites; ++i) {
      for (int j = 0; j < n_sites; ++j) {

        for (auto op_i_str : op_strs) {
          for (auto op_j_str : op_strs) {
            auto op_i = Op(op_i_str, 1.0, i);
            auto op_j = Op(op_j_str, 1.0, j);
            auto op_i_m = matrixC(op_i, block);
            auto op_j_m = matrixC(op_j, block);
            cx_mat anti_comm = op_i_m * op_j_m + op_j_m * op_i_m;

            // Log("{} {}; {} {}", op_i_str, i, op_j_str, j);

            if (((op_i_str == "CDAGUP") && (op_j_str == "CUP")) ||
                ((op_i_str == "CDAGDN") && (op_j_str == "CDN")) ||
                ((op_i_str == "CUP") && (op_j_str == "CDAGUP")) ||
                ((op_i_str == "CDN") && (op_j_str == "CDAGDN"))) {
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

  for (int n_sites = 2; n_sites < 4; ++n_sites) {

    auto block = Electron(n_sites);
    int64_t D = block.size();

    for (int i = 0; i < n_sites; ++i) {

      // Check whether NUMBERUP agrees
      {
        auto op = Op("NUMBERUP", 1.0, i);
        auto m = matrix(op, block);
        auto mi1 = matrix(Op("CDAGUP", 1.0, i), block);
        auto mi2 = matrix(Op("CUP", 1.0, i), block);
        auto mm = mi1 * mi2;
        REQUIRE(norm(mm - m) < 1e-12);
      }

      // Check whether NUMBERDN agrees
      {
        auto op = Op("NUMBERDN", 1.0, i);
        auto m = matrix(op, block);
        auto mi1 = matrix(Op("CDAGDN", 1.0, i), block);
        auto mi2 = matrix(Op("CDN", 1.0, i), block);
        auto mm = mi1 * mi2;
        REQUIRE(norm(mm - m) < 1e-12);
      }

      for (int j = 0; j < n_sites; ++j) {
        if (i != j) {

          // Check whether HOPUP agrees
          {
            auto op = Op("HOPUP", 1.0, {i, j});
            auto m = matrix(op, block);
            auto m1 = matrix(Op("CDAGUP", 1.0, i), block);
            auto m2 = matrix(Op("CUP", 1.0, j), block);
            auto m3 = matrix(Op("CDAGUP", 1.0, j), block);
            auto m4 = matrix(Op("CUP", 1.0, i), block);
            auto mm = -m1 * m2 - m3 * m4;
            REQUIRE(norm(mm - m) < 1e-12);
          }

          // Check whether HOPDN agrees
          {
            auto op = Op("HOPDN", 1.0, {i, j});
            auto m = matrix(op, block);
            auto m1 = matrix(Op("CDAGDN", 1.0, i), block);
            auto m2 = matrix(Op("CDN", 1.0, j), block);
            auto m3 = matrix(Op("CDAGDN", 1.0, j), block);
            auto m4 = matrix(Op("CDN", 1.0, i), block);
            auto mm = -m1 * m2 - m3 * m4;
            REQUIRE(norm(mm - m) < 1e-12);
          }

          // Check whether ISING agrees
          {
            auto op = Op("ISING", 1.0, {i, j});
            auto m = matrix(op, block);

            auto mi1 = matrix(Op("CDAGUP", 1.0, i), block);
            auto mi2 = matrix(Op("CUP", 1.0, i), block);
            auto mi3 = matrix(Op("CDAGDN", 1.0, i), block);
            auto mi4 = matrix(Op("CDN", 1.0, i), block);
            auto szi = 0.5 * mi1 * mi2 - 0.5 * mi3 * mi4;

            auto mj1 = matrix(Op("CDAGUP", 1.0, j), block);
            auto mj2 = matrix(Op("CUP", 1.0, j), block);
            auto mj3 = matrix(Op("CDAGDN", 1.0, j), block);
            auto mj4 = matrix(Op("CDN", 1.0, j), block);
            auto szj = 0.5 * mj1 * mj2 - 0.5 * mj3 * mj4;

            auto mm = szi * szj;
            REQUIRE(norm(mm - m) < 1e-12);
          }

          // Check whether EXCHANGE agrees
          {
            auto op = Op("EXCHANGE", 1.0, {i, j});
            auto m = matrix(op, block);

            auto mi1 = matrix(Op("CDAGUP", 1.0, i), block);
            auto mi2 = matrix(Op("CUP", 1.0, i), block);
            auto mi3 = matrix(Op("CDAGDN", 1.0, i), block);
            auto mi4 = matrix(Op("CDN", 1.0, i), block);
            auto spi = mi1 * mi4;
            auto smi = mi3 * mi2;

            auto mj1 = matrix(Op("CDAGUP", 1.0, j), block);
            auto mj2 = matrix(Op("CUP", 1.0, j), block);
            auto mj3 = matrix(Op("CDAGDN", 1.0, j), block);
            auto mj4 = matrix(Op("CDN", 1.0, j), block);
            auto spj = mj1 * mj4;
            auto smj = mj3 * mj2;

            auto mm = 0.5 * spi * smj + 0.5 * smi * spj;
            REQUIRE(norm(mm - m) < 1e-12);
          }
        }
      }
    }
  }
}
