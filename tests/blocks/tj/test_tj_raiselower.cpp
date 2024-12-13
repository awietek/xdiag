#include "../../catch.hpp"
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/utils/close.hpp>

using namespace xdiag;
using namespace arma;

std::pair<int, int> target_nup_ndn(std::string op_str, int nup, int ndn) {
  if (op_str == "Cdagup") {
    return {nup + 1, ndn};
  } else if (op_str == "Cdagdn") {
    return {nup, ndn + 1};
  } else if (op_str == "Cup") {
    return {nup - 1, ndn};
  } else if (op_str == "Cdn") {
    return {nup, ndn - 1};
  } else {
    Log.err("Invalid op_str");
    return {0, 0};
  }
}

std::pair<int, int> target_nup_ndn(std::string op_str1, std::string op_str2,
                                   int nup, int ndn) {
  int nupd = 0;
  int ndnd = 0;
  if (op_str1 == "Cdagup") {
    ++nupd;
  } else if (op_str1 == "Cdagdn") {
    ++ndnd;
  } else if (op_str1 == "Cup") {
    --nupd;
  } else if (op_str1 == "Cdn") {
    --ndnd;
  } else {
    Log.err("Invalid op_str1");
  }

  if (op_str2 == "Cdagup") {
    ++nupd;
  } else if (op_str2 == "Cdagdn") {
    ++ndnd;
  } else if (op_str2 == "Cup") {
    --nupd;
  } else if (op_str2 == "Cdn") {
    --ndnd;
  } else {
    Log.err("Invalid op_str2");
  }
  return {nup + nupd, ndn + ndnd};
}

bool valid_nup_ndn(int nup, int ndn, int n_sites) {
  return ((nup >= 0) && (ndn >= 0) && (nup + ndn <= n_sites));
}

bool valid_nup_ndn(std::string op_str, int nup, int ndn, int n_sites) {
  auto [nupt, ndnt] = target_nup_ndn(op_str, nup, ndn);
  return valid_nup_ndn(nupt, ndnt, n_sites);
}
bool valid_nup_ndn(std::string op_str1, std::string op_str2, int nup, int ndn,
                   int n_sites) {
  auto [nupt, ndnt] = target_nup_ndn(op_str1, op_str2, nup, ndn);
  return valid_nup_ndn(nupt, ndnt, n_sites);
}

TEST_CASE("tj_raise_lower", "[tj]") try {

  std::vector<std::string> op_strs = {"Cdagup", "Cdagdn", "Cup", "Cdn"};

  for (int n_sites = 1; n_sites < 7; ++n_sites) {
    Log("testing tj anticommutation relations: N={}", n_sites);

    for (int nup = 0; nup <= n_sites; ++nup) {
      for (int ndn = 0; ndn <= n_sites; ++ndn) {
        if (!valid_nup_ndn(nup, ndn, n_sites)) {
          continue;
        }

        auto block = tJ(n_sites, nup, ndn);
        int64_t D = block.size();
        cx_mat id(D, D, fill::eye);
        // Log("n_sites: {}, nup: {}, ndn: {}", n_sites, nup, ndn);

        for (int i = 0; i < n_sites; ++i) {
          for (int j = 0; j < n_sites; ++j) {

            for (auto op_i_str : op_strs) {
              for (auto op_j_str : op_strs) {

                if (!valid_nup_ndn(op_i_str, op_j_str, nup, ndn, n_sites)) {
                  continue;
                }
                if (!valid_nup_ndn(op_i_str, nup, ndn, n_sites)) {
                  continue;
                }
                if (!valid_nup_ndn(op_j_str, nup, ndn, n_sites)) {
                  continue;
                }

                // Log("{} {}; {} {}", op_i_str, i, op_j_str, j);

                auto [nup_i, ndn_i] = target_nup_ndn(op_i_str, nup, ndn);
                auto [nup_j, ndn_j] = target_nup_ndn(op_j_str, nup, ndn);
                auto [nup_ij, ndn_ij] =
                    target_nup_ndn(op_i_str, op_j_str, nup, ndn);

                auto block_i = tJ(n_sites, nup_i, ndn_i);
                auto block_j = tJ(n_sites, nup_j, ndn_j);
                auto block_ij = tJ(n_sites, nup_ij, ndn_ij);

                auto op_i = Op(op_i_str, i);
                auto op_j = Op(op_j_str, j);

                auto op_i_m = matrixC(op_i, block, block_i);
                auto op_ij_m = matrixC(op_j, block_i, block_ij);

                auto op_j_m = matrixC(op_j, block, block_j);
                auto op_ji_m = matrixC(op_i, block_j, block_ij);

                cx_mat anti_comm = op_ji_m * op_j_m + op_ij_m * op_i_m;

                // check the non-fermionic commutation relations of the t-J
                // model
                // see arxiv.org/abs/0706,4236 (tJ model then and now ... Jozef
                // Spalek)
                if (i != j) {
                  REQUIRE(norm(anti_comm) < 1e-12);
                } else {
                  if (((op_i_str == "Cdagup") && (op_j_str == "Cdagup")) ||
                      ((op_i_str == "Cdagup") && (op_j_str == "Cdagdn")) ||
                      ((op_i_str == "Cdagdn") && (op_j_str == "Cdagup")) ||
                      ((op_i_str == "Cdagdn") && (op_j_str == "Cdagdn")) ||
                      ((op_i_str == "Cup") && (op_j_str == "Cup")) ||
                      ((op_i_str == "Cup") && (op_j_str == "Cdn")) ||
                      ((op_i_str == "Cdn") && (op_j_str == "Cup")) ||
                      ((op_i_str == "Cdn") && (op_j_str == "Cdn"))) {
                    REQUIRE(norm(anti_comm) < 1e-12);
                  } else if (((op_i_str == "Cdagup") && (op_j_str == "Cup")) ||
                             ((op_i_str == "Cup") && (op_j_str == "Cdagup"))) {
                    auto ndn_op = Op("Ndn", i);
                    auto ndn_op_m = matrix(ndn_op, block);
                    REQUIRE(norm(anti_comm + ndn_op_m - id) < 1e-12);
                  } else if (((op_i_str == "Cdagdn") && (op_j_str == "Cdn")) ||
                             ((op_i_str == "Cdn") && (op_j_str == "Cdagdn"))) {
                    auto nup_op = Op("Nup", i);
                    auto nup_op_m = matrix(nup_op, block);
                    REQUIRE(norm(anti_comm + nup_op_m - id) < 1e-12);
                  } else if ((op_i_str == "Cdagup") && (op_j_str == "Cdn")) {
                    REQUIRE(norm(anti_comm - op_ji_m * op_j_m) < 1e-12);
                  } else if ((op_i_str == "Cdn") && (op_j_str == "Cdagup")) {
                    REQUIRE(norm(anti_comm - op_ij_m * op_i_m) < 1e-12);
                  } else if ((op_i_str == "Cdagdn") && (op_j_str == "Cup")) {
                    REQUIRE(norm(anti_comm - op_ji_m * op_j_m) < 1e-12);
                  } else if ((op_i_str == "Cup") && (op_j_str == "Cdagdn")) {
                    REQUIRE(norm(anti_comm - op_ij_m * op_i_m) < 1e-12);
                  } else {
                    Log("unchecked"); // shouldn't happen
                  }
                }

              } // loop op_i_str
            } // loop op_j_str

          } // loop i
        } // loop j

      } // loop nup
    } // loop ndn

  } // loop n_sites

  // for (int n_sites = 2; n_sites < 3; ++n_sites) {

  //   auto block = tJ(n_sites);
  //   int64_t D = block.size();

  //   for (int i = 0; i < n_sites; ++i) {

  //     // Check whether Nup agrees
  //     {
  //       auto op = Op("Nup", i);
  //       auto m = matrix(op, block);
  //       auto mi1 = matrix(Op("Cdagup", i), block);
  //       auto mi2 = matrix(Op("Cup", i), block);
  //       auto mm = mi1 * mi2;
  //       REQUIRE(norm(mm - m) < 1e-12);
  //     }

  //     // Check whether Ndn agrees
  //     {
  //       auto op = Op("Ndn", i);
  //       auto m = matrix(op, block);
  //       auto mi1 = matrix(Op("Cdagdn", i), block);
  //       auto mi2 = matrix(Op("Cdn", i), block);
  //       auto mm = mi1 * mi2;
  //       REQUIRE(norm(mm - m) < 1e-12);
  //     }

  //     for (int j = 0; j < n_sites; ++j) {
  //       if (i != j) {

  //         // Check whether Hopup agrees
  //         {
  //           auto op = Op("Hopup", {i, j});
  //           auto m = matrix(op, block);
  //           auto m1 = matrix(Op("Cdagup", i), block);
  //           auto m2 = matrix(Op("Cup", j), block);
  //           auto m3 = matrix(Op("Cdagup", j), block);
  //           auto m4 = matrix(Op("Cup", i), block);
  //           auto mm = -m1 * m2 - m3 * m4;
  //           REQUIRE(norm(mm - m) < 1e-12);
  //         }

  //         // Check whether Hopdn agrees
  //         {
  //           auto op = Op("Hopdn", {i, j});
  //           auto m = matrix(op, block);
  //           auto m1 = matrix(Op("Cdagdn", i), block);
  //           auto m2 = matrix(Op("Cdn", j), block);
  //           auto m3 = matrix(Op("Cdagdn", j), block);
  //           auto m4 = matrix(Op("Cdn", i), block);
  //           auto mm = -m1 * m2 - m3 * m4;
  //           REQUIRE(norm(mm - m) < 1e-12);
  //         }

  //         // Check whether SzSz agrees
  //         {
  //           auto op = Op("SzSz", {i, j});
  //           auto m = matrix(op, block);

  //           auto mi1 = matrix(Op("Cdagup", i), block);
  //           auto mi2 = matrix(Op("Cup", i), block);
  //           auto mi3 = matrix(Op("Cdagdn", i), block);
  //           auto mi4 = matrix(Op("Cdn", i), block);
  //           auto szi = 0.5 * mi1 * mi2 - 0.5 * mi3 * mi4;

  //           auto mj1 = matrix(Op("Cdagup", j), block);
  //           auto mj2 = matrix(Op("Cup", j), block);
  //           auto mj3 = matrix(Op("Cdagdn", j), block);
  //           auto mj4 = matrix(Op("Cdn", j), block);
  //           auto szj = 0.5 * mj1 * mj2 - 0.5 * mj3 * mj4;

  //           auto mm = szi * szj;
  //           REQUIRE(norm(mm - m) < 1e-12);
  //         }

  //         // Check whether Exchange agrees
  //         {
  //           auto op = Op("Exchange", {i, j});
  //           auto m = matrix(op, block);

  //           auto mi1 = matrix(Op("Cdagup", i), block);
  //           auto mi2 = matrix(Op("Cup", i), block);
  //           auto mi3 = matrix(Op("Cdagdn", i), block);
  //           auto mi4 = matrix(Op("Cdn", i), block);
  //           auto spi = mi1 * mi4;
  //           auto smi = mi3 * mi2;

  //           auto mj1 = matrix(Op("Cdagup", j), block);
  //           auto mj2 = matrix(Op("Cup", j), block);
  //           auto mj3 = matrix(Op("Cdagdn", j), block);
  //           auto mj4 = matrix(Op("Cdn", j), block);
  //           auto spj = mj1 * mj4;
  //           auto smj = mj3 * mj2;

  //           auto mm = 0.5 * spi * smj + 0.5 * smi * spj;
  //           REQUIRE(norm(mm - m) < 1e-12);
  //         }
  //       }
  //     }
  //   }
  // }
} catch (xdiag::Error const &e) {
  error_trace(e);
}
