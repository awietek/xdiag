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

static OpSum bose_hubbard_opsum(int64_t nsites, double t, double U) {

  auto ops = OpSum();
  for (int i = 0; i < nsites; ++i) {
    ops += t * Op("Hop", {i, (i + 1) % nsites});
  }
  ops += U * Op("HubbardU");
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
    int64_t d = 5;
    auto block = Boson(nsites, d);

    Log("Bose-Hubbard tests");
    {
      auto [e0, psi0] = eig0(bose_hubbard_opsum(nsites, 1.0, 1.0), block);
      double e0_dmrg = -9.354167806978255;
      Log("e0: {:.16f}, e0_dmrg: {:.16f}", e0, e0_dmrg);
      REQUIRE(isapprox(e0, e0_dmrg, 1e-8, 1e-8));

      arma::mat AdagA_dmrg = {
          {2.0, 1.87457, 1.82526, 1.87457},
          {1.87457, 2.0, 1.87457, 1.82526},
          {1.82526, 1.87457, 2.0, 1.87457},
          {1.87457, 1.82526, 1.87457, 2.0},
      };
      auto AdagA = correlation_matrix(psi0, "Adag", "A");
      REQUIRE(isapprox(AdagA_dmrg, AdagA, 1e-5, 1e-5));
    }
    {
      auto [e0, psi0] = eig0(bose_hubbard_opsum(nsites, 1.0, 2.0), block);
      double e0_dmrg = -5.921611286247851;
      Log("e0: {:.16f}, e0_dmrg: {:.16f}", e0, e0_dmrg);
      REQUIRE(isapprox(e0, e0_dmrg, 1e-8, 1e-8));

      arma::mat AdagA_dmrg = {{1.25, 1.17416, 1.14403, 1.17416},
                              {1.17416, 1.25, 1.17416, 1.14403},
                              {1.14403, 1.17416, 1.25, 1.17416},
                              {1.17416, 1.14403, 1.17416, 1.25}};
      auto AdagA = correlation_matrix(psi0, "Adag", "A");
      REQUIRE(isapprox(AdagA_dmrg, AdagA, 1e-5, 1e-5));
    }

    {
      auto [e0, psi0] = eig0(bose_hubbard_opsum(nsites, 1.0, 4.0), block);
      double e0_dmrg = -4.116103420851046;
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
      auto [e0, psi0] = eig0(bose_hubbard_opsum(nsites, 1.0, 8.0), block);
      double e0_dmrg = -3.32638322925368;
      Log("e0: {:.16f}, e0_dmrg: {:.16f}", e0, e0_dmrg);
      REQUIRE(isapprox(e0, e0_dmrg, 1e-8, 1e-8));

      arma::mat AdagA_dmrg = {{0.75, 0.547103, 0.501433, 0.547103},
                              {0.547103, 0.75, 0.547103, 0.501433},
                              {0.501433, 0.547103, 0.75, 0.547103},
                              {0.547103, 0.501433, 0.547103, 0.75}};
      auto AdagA = correlation_matrix(psi0, "Adag", "A");
      REQUIRE(isapprox(AdagA_dmrg, AdagA, 1e-5, 1e-5));
    }
  }

  // Test whether quantum numbers work
  for (int nsites = 2; nsites < 5; ++nsites) {
    Log("Bose-Hubbard chain symmetry test: N = {}", nsites);
    auto ops = bose_hubbard_opsum(nsites, 1.0, 2.0);
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
