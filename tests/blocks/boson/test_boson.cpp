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
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

TEST_CASE("boson", "[boson]") try {
  using namespace xdiag;
  int64_t nsites = 4;
  int64_t d = 5;
  auto block = Boson(nsites, d);
  // for (auto state : block) {
  //   Log(to_string(state, block));
  // }
  Log("Bose-Hubbard tests");
  
  double t = 1.0;
  auto ops = OpSum();
  for (int i = 0; i < nsites; ++i) {
    ops += t * Op("Hop", {i, (i + 1) % nsites});
  }
  ops += "U" * Op("HubbardU");

  {
    ops["U"] = 1.0;
    auto [e0, psi0] = eig0(ops, block);
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
    ops["U"] = 2.0;
    auto [e0, psi0] = eig0(ops, block);
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
    ops["U"] = 4.0;
    auto [e0, psi0] = eig0(ops, block);
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
    ops["U"] = 8.0;
    auto [e0, psi0] = eig0(ops, block);
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
} catch (xdiag::Error e) {
  error_trace(e);
}
