// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <vector>

#include <tests/catch.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

// Validates the electron matrix kernels (and their Jordan-Wigner signs) by
// reconstructing each named operator from products of the elementary
// Cdagup/Cup/Cdagdn/Cdn matrices on the full (number-non-conserving) block.
TEST_CASE("electron_raise_lower", "[electron]") try {
  using namespace arma;

  // Fermionic anti-commutation relations on the full Fock space.
  std::vector<std::string> op_strs = {"Cdagup", "Cdagdn", "Cup", "Cdn"};
  for (int nsites = 1; nsites < 4; ++nsites) {
    Log("electron anti-commutation relations (full matrix): N={}", nsites);
    auto block = Electron(nsites);
    int64_t D = block.size();
    cx_mat id(D, D, fill::eye);
    cx_mat zeros(D, D, fill::zeros);

    for (int i = 0; i < nsites; ++i) {
      for (int j = 0; j < nsites; ++j) {
        for (auto op_i_str : op_strs) {
          for (auto op_j_str : op_strs) {
            cx_mat op_i_m = matrixC(Op(op_i_str, i), block);
            cx_mat op_j_m = matrixC(Op(op_j_str, j), block);
            cx_mat anti_comm = op_i_m * op_j_m + op_j_m * op_i_m;

            bool conjugate_pair =
                ((op_i_str == "Cdagup") && (op_j_str == "Cup")) ||
                ((op_i_str == "Cdagdn") && (op_j_str == "Cdn")) ||
                ((op_i_str == "Cup") && (op_j_str == "Cdagup")) ||
                ((op_i_str == "Cdn") && (op_j_str == "Cdagdn"));
            if (conjugate_pair && (i == j)) {
              REQUIRE(isapprox(anti_comm, id));
            } else {
              REQUIRE(isapprox(anti_comm, zeros));
            }
          }
        }
      }
    }
  }

  // Named operators reconstructed from elementary operators (i != j).
  for (int nsites = 2; nsites < 4; ++nsites) {
    Log("electron named operators vs elementary products: N={}", nsites);
    auto block = Electron(nsites);

    for (int i = 0; i < nsites; ++i) {

      // Nup{i} = Cdagup{i} Cup{i}
      {
        mat m = matrix(Op("Nup", i), block);
        mat mm = matrix(Op("Cdagup", i), block) * matrix(Op("Cup", i), block);
        REQUIRE(isapprox(mm, m));
      }
      // Ndn{i} = Cdagdn{i} Cdn{i}
      {
        mat m = matrix(Op("Ndn", i), block);
        mat mm = matrix(Op("Cdagdn", i), block) * matrix(Op("Cdn", i), block);
        REQUIRE(isapprox(mm, m));
      }

      for (int j = 0; j < nsites; ++j) {
        if (i == j) {
          continue;
        }

        // Hopup{i,j} = -Cdagup{i} Cup{j} - Cdagup{j} Cup{i}
        {
          mat m = matrix(Op("Hopup", {i, j}), block);
          mat m1 = matrix(Op("Cdagup", i), block);
          mat m2 = matrix(Op("Cup", j), block);
          mat m3 = matrix(Op("Cdagup", j), block);
          mat m4 = matrix(Op("Cup", i), block);
          REQUIRE(isapprox(mat(-m1 * m2 - m3 * m4), m));
        }

        // Hopdn{i,j} = -Cdagdn{i} Cdn{j} - Cdagdn{j} Cdn{i}
        {
          mat m = matrix(Op("Hopdn", {i, j}), block);
          mat m1 = matrix(Op("Cdagdn", i), block);
          mat m2 = matrix(Op("Cdn", j), block);
          mat m3 = matrix(Op("Cdagdn", j), block);
          mat m4 = matrix(Op("Cdn", i), block);
          REQUIRE(isapprox(mat(-m1 * m2 - m3 * m4), m));
        }

        // SzSz{i,j} = Sz_i Sz_j, Sz_s = 1/2 (Cdagup_s Cup_s - Cdagdn_s Cdn_s)
        {
          mat m = matrix(Op("SzSz", {i, j}), block);
          mat szi = 0.5 * matrix(Op("Cdagup", i), block) *
                        matrix(Op("Cup", i), block) -
                    0.5 * matrix(Op("Cdagdn", i), block) *
                        matrix(Op("Cdn", i), block);
          mat szj = 0.5 * matrix(Op("Cdagup", j), block) *
                        matrix(Op("Cup", j), block) -
                    0.5 * matrix(Op("Cdagdn", j), block) *
                        matrix(Op("Cdn", j), block);
          REQUIRE(isapprox(mat(szi * szj), m));
        }

        // Exchange{i,j} = 1/2 (S+_i S-_j + S-_i S+_j), with
        // S+_s = Cdagup_s Cdn_s, S-_s = Cdagdn_s Cup_s
        {
          mat m = matrix(Op("Exchange", {i, j}), block);
          mat spi = matrix(Op("Cdagup", i), block) * matrix(Op("Cdn", i), block);
          mat smi = matrix(Op("Cdagdn", i), block) * matrix(Op("Cup", i), block);
          mat spj = matrix(Op("Cdagup", j), block) * matrix(Op("Cdn", j), block);
          mat smj = matrix(Op("Cdagdn", j), block) * matrix(Op("Cup", j), block);
          REQUIRE(isapprox(mat(0.5 * spi * smj + 0.5 * smi * spj), m));
        }

        // NupNdn / NdnNup / NupNup / NdnNdn
        {
          mat nupi = matrix(Op("Nup", i), block);
          mat ndni = matrix(Op("Ndn", i), block);
          mat nupj = matrix(Op("Nup", j), block);
          mat ndnj = matrix(Op("Ndn", j), block);
          REQUIRE(isapprox(mat(nupi * ndnj), matrix(Op("NupNdn", {i, j}), block)));
          REQUIRE(isapprox(mat(ndni * nupj), matrix(Op("NdnNup", {i, j}), block)));
          REQUIRE(isapprox(mat(nupi * nupj), matrix(Op("NupNup", {i, j}), block)));
          REQUIRE(isapprox(mat(ndni * ndnj), matrix(Op("NdnNdn", {i, j}), block)));
        }
      }
    }
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// Spin operators S+, S-, Sx, Sy (which expand to elementary-operator strings)
// and the canonical su(2) commutation relations, on the full Fock space:
//   S+_s = Cdagup_s Cdn_s,  S-_s = Cdagdn_s Cup_s,  Sz_s = 1/2 (Nup_s - Ndn_s)
//   Sx_s = 1/2 (S+_s + S-_s),  Sy_s = -i/2 (S+_s - S-_s)
//   [S+_i, S-_j] = 2 Sz_i delta_ij,  [Sz_i, S+-_j] = +- S+-_i delta_ij
TEST_CASE("electron_commutators", "[electron]") try {
  using namespace arma;
  complex I(0.0, 1.0);

  for (int nsites = 2; nsites < 4; ++nsites) {
    Log("electron spin operators and su(2) commutators: N={}", nsites);
    auto block = Electron(nsites);
    int64_t D = block.size();
    mat zeros(D, D, fill::zeros);

    // S+, S-, Sx, Sy reconstructed from the elementary operators.
    for (int s = 0; s < nsites; ++s) {
      mat sp_ref = matrix(Op("Cdagup", s), block) * matrix(Op("Cdn", s), block);
      mat sm_ref = matrix(Op("Cdagdn", s), block) * matrix(Op("Cup", s), block);
      REQUIRE(isapprox(matrix(Op("S+", s), block), sp_ref));
      REQUIRE(isapprox(matrix(Op("S-", s), block), sm_ref));
      REQUIRE(isapprox(matrix(Op("Sx", s), block), mat(0.5 * (sp_ref + sm_ref))));

      cx_mat spc =
          matrixC(Op("Cdagup", s), block) * matrixC(Op("Cdn", s), block);
      cx_mat smc =
          matrixC(Op("Cdagdn", s), block) * matrixC(Op("Cup", s), block);
      REQUIRE(isapprox(matrixC(Op("Sy", s), block),
                       cx_mat(-0.5 * I * spc + 0.5 * I * smc)));
    }

    // Canonical su(2) commutation relations (using the named operators).
    for (int i = 0; i < nsites; ++i) {
      mat sz_i = matrix(Op("Sz", i), block);
      for (int j = 0; j < nsites; ++j) {
        mat sp_i = matrix(Op("S+", i), block);
        mat sm_i = matrix(Op("S-", i), block);
        mat sp_j = matrix(Op("S+", j), block);
        mat sm_j = matrix(Op("S-", j), block);

        mat comm_pm = sp_i * sm_j - sm_j * sp_i; // [S+_i, S-_j]
        mat comm_zp = sz_i * sp_j - sp_j * sz_i; // [Sz_i, S+_j]
        mat comm_zm = sz_i * sm_j - sm_j * sz_i; // [Sz_i, S-_j]

        if (i == j) {
          REQUIRE(isapprox(comm_pm, mat(2.0 * sz_i)));
          REQUIRE(isapprox(comm_zp, sp_i));
          REQUIRE(isapprox(comm_zm, mat(-sm_i)));
        } else {
          REQUIRE(isapprox(comm_pm, zeros));
          REQUIRE(isapprox(comm_zp, zeros));
          REQUIRE(isapprox(comm_zm, zeros));
        }
      }
    }
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
