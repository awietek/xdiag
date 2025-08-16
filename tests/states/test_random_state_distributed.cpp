// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include "../blocks/electron/testcases_electron.hpp"
#include <xdiag/algebra/algebra.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/utils/xdiag_show.hpp>

#include <iostream>
#include <set>

using namespace xdiag;

TEST_CASE("random_state", "[states]") try {
  using namespace xdiag::testcases::electron;
  {
    using namespace arma;

    Log("random matrix state");

    auto b = SpinhalfDistributed(6, 3);
    int64_t dim = size(b);
    int64_t ncols = 3;
    int64_t seed = 42;

    auto rstate = random_state(b, true, ncols, seed, true);
    auto ovlp = matrix_dot(rstate, rstate);
    REQUIRE(approx_equal(ovlp, mat(eye<mat>(ncols, ncols)), "absdiff", 1e-14));

    auto rstateC = random_state(b, false, ncols, seed, true);
    auto ovlpC = matrix_dotC(rstateC, rstateC);
    REQUIRE(approx_equal(ovlpC, cx_mat(eye<cx_mat>(ncols, ncols)), "absdiff",
                         1e-14));
  }
}
XDIAG_CATCH
