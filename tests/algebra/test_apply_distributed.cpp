// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/operators/logic/order.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/operators/logic/symmetrize.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>

using namespace xdiag;
using namespace arma;

void test_apply_matrix_distributed(OpSum ops, Block block, bool real, int64_t ncols) try {
  int64_t seed = 42;
  bool normalized = true;

  auto r = random_state(block, real, ncols, seed, normalized);

  auto Ar = apply(ops, r);

  auto r2 = random_state(Ar.block(), real, ncols);
  

  if (real && isreal(ops)) {
    auto rAr = matrix_dot(r2, Ar);
    for (int i = 0; i < ncols; ++i) {
      for (int j = 0; j < ncols; ++j) {
        double rAr_ij = dot(r2.col(i), apply(ops, r.col(j)));
        REQUIRE(isapprox(rAr_ij, rAr(i, j)));
      }
    }
  } else {
    auto rAr = matrix_dotC(r2, Ar);
    for (int i = 0; i < ncols; ++i) {
      for (int j = 0; j < ncols; ++j) {
        complex rAr_ij = dotC(r2.col(i), apply(ops, r.col(j)));
        REQUIRE(isapprox(rAr_ij, rAr(i, j)));
      }
    }
  }
}
XDIAG_CATCH

TEST_CASE("apply_distributed", "[algebra]") try {

  {
    OpSum ops;
    ops += Op("Hop", {0, 1});
    ElectronDistributed block(4, 2, 2);
    test_apply_matrix_distributed(ops, block, true, 3);
    test_apply_matrix_distributed(ops, block, true, 4);
    test_apply_matrix_distributed(ops, block, false, 3);
    test_apply_matrix_distributed(ops, block, false, 4);
  }

  {
    OpSum ops;
    ops += Op("Cdagup", 1);
    ElectronDistributed block(4, 2, 2);
    test_apply_matrix_distributed(ops, block, true, 3);
    test_apply_matrix_distributed(ops, block, true, 4);
    test_apply_matrix_distributed(ops, block, false, 3);
    test_apply_matrix_distributed(ops, block, false, 4);
  }
} catch (Error const &e) {
  error_trace(e);
}

