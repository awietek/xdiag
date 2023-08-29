#include "../../catch.hpp"

#include <hydra/blocks/electron/electron_matrix.h>
#include <hydra/blocks/electron/electron_apply.h>
#include <hydra/algebra/algebra.h>
#include <hydra/algebra/matrix.h>
#include <hydra/algorithms/sparse_diag.h>
#include <hydra/utils/close.h>

TEST_CASE("electron_raise_lower", "[electron]") {
  using namespace hydra;
  using namespace arma;
  std::vector<std::string> op_strs = {"CDAGUP", "CDAGDN", "CUP", "CDN"};

  for (int n_sites = 1; n_sites < 4; ++n_sites) {
    Log("testing electron anticommutation relations (full matrix): N={}",
        n_sites);

    auto block = Electron(n_sites);
    idx_t D = block.size();

    cx_mat id(D, D, fill::eye);

    for (int i = 0; i < n_sites; ++i) {
      for (int j = 0; j < n_sites; ++j) {

        for (auto op_i_str : op_strs) {
          for (auto op_j_str : op_strs) {
            auto op_i = Bond(op_i_str, i);
            auto op_j = Bond(op_j_str, j);
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
    idx_t D = block.size();

    for (int i = 0; i < n_sites; ++i) {

      // Check whether NUMBERUP agrees
      {
        auto bond = Bond("NUMBERUP", i);
        auto m = matrix(bond, block);
        auto mi1 = matrix(Bond("CDAGUP", i), block);
        auto mi2 = matrix(Bond("CUP", i), block);
        auto mm = mi1 * mi2;
        REQUIRE(norm(mm - m) < 1e-12);
      }

      // Check whether NUMBERDN agrees
      {
        auto bond = Bond("NUMBERDN", i);
        auto m = matrix(bond, block);
        auto mi1 = matrix(Bond("CDAGDN", i), block);
        auto mi2 = matrix(Bond("CDN", i), block);
        auto mm = mi1 * mi2;
        REQUIRE(norm(mm - m) < 1e-12);
      }

      for (int j = 0; j < n_sites; ++j) {
        if (i != j) {

          // Check whether HOPUP agrees
          {
            auto bond = Bond("HOPUP", {i, j});
            auto m = matrix(bond, block);
            auto m1 = matrix(Bond("CDAGUP", i), block);
            auto m2 = matrix(Bond("CUP", j), block);
            auto m3 = matrix(Bond("CDAGUP", j), block);
            auto m4 = matrix(Bond("CUP", i), block);
            auto mm = -m1 * m2 - m3 * m4;
            REQUIRE(norm(mm - m) < 1e-12);
          }

          // Check whether HOPDN agrees
          {
            auto bond = Bond("HOPDN", {i, j});
            auto m = matrix(bond, block);
            auto m1 = matrix(Bond("CDAGDN", i), block);
            auto m2 = matrix(Bond("CDN", j), block);
            auto m3 = matrix(Bond("CDAGDN", j), block);
            auto m4 = matrix(Bond("CDN", i), block);
            auto mm = -m1 * m2 - m3 * m4;
            REQUIRE(norm(mm - m) < 1e-12);
          }

          // Check whether ISING agrees
          {
            auto bond = Bond("ISING", {i, j});
            auto m = matrix(bond, block);

            auto mi1 = matrix(Bond("CDAGUP", i), block);
            auto mi2 = matrix(Bond("CUP", i), block);
            auto mi3 = matrix(Bond("CDAGDN", i), block);
            auto mi4 = matrix(Bond("CDN", i), block);
            auto szi = 0.5 * mi1 * mi2 - 0.5 * mi3 * mi4;

            auto mj1 = matrix(Bond("CDAGUP", j), block);
            auto mj2 = matrix(Bond("CUP", j), block);
            auto mj3 = matrix(Bond("CDAGDN", j), block);
            auto mj4 = matrix(Bond("CDN", j), block);
            auto szj = 0.5 * mj1 * mj2 - 0.5 * mj3 * mj4;

            auto mm = szi * szj;
            REQUIRE(norm(mm - m) < 1e-12);
          }

          // Check whether EXCHANGE agrees
          {
            auto bond = Bond("EXCHANGE", {i, j});
            auto m = matrix(bond, block);

            auto mi1 = matrix(Bond("CDAGUP", i), block);
            auto mi2 = matrix(Bond("CUP", i), block);
            auto mi3 = matrix(Bond("CDAGDN", i), block);
            auto mi4 = matrix(Bond("CDN", i), block);
            auto spi = mi1 * mi4;
            auto smi = mi3 * mi2;

            auto mj1 = matrix(Bond("CDAGUP", j), block);
            auto mj2 = matrix(Bond("CUP", j), block);
            auto mj3 = matrix(Bond("CDAGDN", j), block);
            auto mj4 = matrix(Bond("CDN", j), block);
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
