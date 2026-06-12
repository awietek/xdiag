// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/matrices/matrix.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/correlation_matrix.hpp>
#include <xdiag/states/expect.hpp>
#include <xdiag/states/inner.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

using namespace xdiag;

static OpSum bose_hubbard_opsum(int64_t nsites, double t, double U, double mu) {

  auto ops = OpSum();
  for (int i = 0; i < nsites; ++i) {
    ops += t * Op("Hop", {i, (i + 1) % nsites});
  }
  ops += U * Op("HubbardU");
  ops += mu * Op("TotalN");
  return ops;
}

TEST_CASE("boson", "[boson]") try {
  {
    Boson block;
    REQUIRE_THROWS(block = Boson(3, 3, 9));
    REQUIRE_THROWS(block = Boson(-3, 3, 9));
    REQUIRE_THROWS(block = Boson(-3, 1, 9));
    REQUIRE_THROWS(block = Boson(3, 3, -9));
    REQUIRE_THROWS(block = Boson(3, 3333, 2));
  }
  {
    int64_t nsites = 4;
    double t = 1.0;
    double mu = 0.5;
    int64_t d = 5;
    auto block = Boson(nsites, d);

    Log("Bose-Hubbard tests");
    {
      double U = 1.0;
      auto [e0, psi0] = eig0(bose_hubbard_opsum(nsites, t, U, mu), block);
      double e0_dmrg = -5.660431125405925;
      Log("e0: {:.16f}, e0_dmrg: {:.16f}", e0, e0_dmrg);
      REQUIRE(isapprox(e0, e0_dmrg, 1e-8, 1e-8));

      arma::mat AdagA_dmrg = {
          {1.75, 1.67525, 1.6428, 1.67525},
          {1.67525, 1.75, 1.67525, 1.6428},
          {1.6428, 1.67525, 1.75, 1.67525},
          {1.67525, 1.6428, 1.67525, 1.75},
      };
      auto AdagA = correlation_matrix(psi0, "Adag", "A");
      REQUIRE(isapprox(AdagA_dmrg, AdagA, 1e-5, 1e-5));
    }
    {
      double U = 2.0;
      auto [e0, psi0] = eig0(bose_hubbard_opsum(nsites, t, U, mu), block);
      double e0_dmrg = -3.6138649391176814;
      Log("e0: {:.16f}, e0_dmrg: {:.16f}", e0, e0_dmrg);
      REQUIRE(isapprox(e0, e0_dmrg, 1e-8, 1e-8));

      arma::mat AdagA_dmrg = {{1.0, 0.945091, 0.921613, 0.945091},
                              {0.945091, 1.0, 0.945091, 0.921613},
                              {0.921613, 0.945091, 1.0, 0.945091},
                              {0.945091, 0.921613, 0.945091, 1.0}};
      auto AdagA = correlation_matrix(psi0, "Adag", "A");
      REQUIRE(isapprox(AdagA_dmrg, AdagA, 1e-5, 1e-5));
    }

    {
      double U = 4.0;
      auto [e0, psi0] = eig0(bose_hubbard_opsum(nsites, t, U, mu), block);
      double e0_dmrg = -2.6161034208510463;
      Log("e0: {:.16f}, e0_dmrg: {:.16f}", e0, e0_dmrg);
      REQUIRE(isapprox(e0, e0_dmrg, 1e-8, 1e-8));

      arma::mat AdagA_dmrg = {{0.75, 0.659022, 0.625756, 0.659022},
                              {0.659022, 0.75, 0.659022, 0.625756},
                              {0.625756, 0.659022, 0.75, 0.659022},
                              {0.659022, 0.625756, 0.659022, 0.75}};
      auto AdagA = correlation_matrix(psi0, "Adag", "A");
      REQUIRE(isapprox(AdagA_dmrg, AdagA, 1e-5, 1e-5));
    }

    {
      double U = 8.0;
      auto [e0, psi0] = eig0(bose_hubbard_opsum(nsites, t, U, mu), block);
      double e0_dmrg = -2.207750943219353;
      Log("e0: {:.16f}, e0_dmrg: {:.16f}", e0, e0_dmrg);
      REQUIRE(isapprox(e0, e0_dmrg, 1e-8, 1e-8));

      arma::mat AdagA_dmrg = {{0.5, 0.43556, 0.392426, 0.43556},
                              {0.43556, 0.5, 0.43556, 0.392426},
                              {0.392426, 0.43556, 0.5, 0.43556},
                              {0.43556, 0.392426, 0.43556, 0.5}};
      auto AdagA = correlation_matrix(psi0, "Adag", "A");
      REQUIRE(isapprox(AdagA_dmrg, AdagA, 1e-5, 1e-5));
    }
  }

  // Test whether quantum numbers work
  for (int nsites = 2; nsites < 5; ++nsites) {
    Log("Bose-Hubbard chain symmetry test: N = {}", nsites);
    auto ops = bose_hubbard_opsum(nsites, 1.0, 2.0, 0.0);
    for (int d = 2; d < 6; ++d) {
      auto block = Boson(nsites, d);
      double e0 = eigval0(ops, block);
      // Log("e0: {:.6f}", e0);
      double e0nmin = 999.999999;

      for (int number = 0; number < nsites * (d - 1); ++number) {
        auto blockn = Boson(nsites, d, number);
        auto [e0n, psi0n] = eig0(ops, blockn);
        // Log("n: {}, e0n: {:.6f}", number, e0n);
        if (e0n < e0nmin) {
          e0nmin = e0n;
        }

        auto AdagAn = correlation_matrixC(psi0n, "Adag", "A");

        // XDIAG_SHOW(AdagAn);

        double e0nkmin = 99999.999;
        int kmin = 0;
        for (int k = 0; k < nsites; ++k) {
          auto blocknk =
              Boson(nsites, d, number, cyclic_group_irrep(nsites, k));
          // XDIAG_SHOW(blocknk);
          if (dim(blocknk) > 0) {
            auto [e0nk, psi0nk] = eig0(ops, blocknk);
            // Log("n: {}, k: {}, e0nk: {:.6f}", number, k, e0nk);
            if (e0nk < e0nkmin) {
              e0nkmin = e0nk;
              kmin = k;
            }
          }
        }
        auto blocknk =
            Boson(nsites, d, number, cyclic_group_irrep(nsites, kmin));
        auto [e0nk, psi0nk] = eig0(ops, blocknk);
        auto AdagAnk = correlation_matrixC(psi0nk, "Adag", "A");
        // XDIAG_SHOW(AdagAnk);

        REQUIRE(isapprox(e0n, e0nkmin));
        REQUIRE(isapprox(AdagAn, AdagAnk, 1e-6, 1e-6));
      }
      REQUIRE(isapprox(e0, e0nmin));

      double e0kmin = 99999.999;
      for (int k = 0; k < nsites; ++k) {
        auto blockk = Boson(nsites, d, cyclic_group_irrep(nsites, k));
        if (dim(blockk) > 0) {
          double e0k = eigval0(ops, blockk);
          // Log("k: {}, e0k: {:.6f}", k, e0k);
          if (e0k < e0kmin) {
            e0kmin = e0k;
          }
        }
      }
      REQUIRE(isapprox(e0, e0kmin));
      // Log("");
    }
  }

} catch (xdiag::Error e) {
  error_trace(e);
}

TEST_CASE("bosoncommutation", "[boson]") try {
  Log("Boson (truncated) commutation relations");

  // Bosons on a Fock space truncated to occupations 0..d-1 satisfy
  //   [a_i, a_j]   = 0
  //   [adag_i, adag_j] = 0
  //   [a_i, adag_j]    = 0                         (i != j)
  //   [a_i, adag_i]    = I - d P_{n=d-1}           (truncated!)
  // The truncation modifies the canonical [a, adag] = I only on the highest
  // occupation n = d-1, where adag|d-1> = 0, giving the diagonal entry -(d-1)
  // instead of +1.
  for (int nsites = 2; nsites <= 3; ++nsites) {
    for (int d = 2; d <= 4; ++d) {
      auto block = Boson(nsites, d); // full Fock space (square matrices)
      int64_t D = block.dim();
      arma::mat zeros(D, D, arma::fill::zeros);

      for (int i = 0; i < nsites; ++i) {
        for (int j = 0; j < nsites; ++j) {
          arma::mat adagi = matrix(Op("Adag", i), block);
          arma::mat adagj = matrix(Op("Adag", j), block);
          arma::mat ai = matrix(Op("A", i), block);
          arma::mat aj = matrix(Op("A", j), block);

          arma::mat aa = ai * aj - aj * ai; // [a_i, a_j]
          arma::mat adagadag =
              adagi * adagj - adagj * adagi;          // [adag_i, adag_j]
          arma::mat a_adag = ai * adagj - adagj * ai; // [a_i, adag_j]
          arma::mat adag_a = adagi * aj - aj * adagi; // [adag_i, a_j]

          REQUIRE(isapprox(aa, zeros));
          REQUIRE(isapprox(adagadag, zeros));

          if (i != j) {
            REQUIRE(isapprox(a_adag, zeros));
            REQUIRE(isapprox(adag_a, zeros));
          } else {
            // Truncated commutator built from the number operator: the diagonal
            // is +1 except -(d-1) on states with maximal occupation at site i.
            arma::mat ni = matrix(Op("N", i), block); // diagonal n_i
            arma::vec nd = ni.diag();
            arma::mat expected(D, D, arma::fill::zeros);
            for (int64_t k = 0; k < D; ++k) {
              bool maxed = nd(k) > static_cast<double>(d) - 1.5; // n_i == d-1
              expected(k, k) = maxed ? -static_cast<double>(d - 1) : 1.0;
            }
            REQUIRE(isapprox(a_adag, expected));
            REQUIRE(isapprox(adag_a, arma::mat(-expected)));
          }
        }
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
}
