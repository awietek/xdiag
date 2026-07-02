// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <vector>

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/operators/valid.hpp>
#include <xdiag/utils/error.hpp>

// Negative tests for the Op validation layer. Every branch here is a rejection
// path that is only reachable with malformed input, so the positive-path tests
// never touch them.
TEST_CASE("check_valid Op", "[operators]") try {
  using namespace xdiag;
  using namespace xdiag::operators;

  // --- Valid Ops must not throw ---
  {
    REQUIRE_NOTHROW(check_valid(Op("Id")));            // site-free, no matrix
    REQUIRE_NOTHROW(check_valid(Op("Sz", int64_t(0)))); // single site
    REQUIRE_NOTHROW(
        check_valid(Op("SdotS", std::vector<int64_t>{0, 1}))); // two sites
    REQUIRE_NOTHROW(check_valid(Op("HubbardU"))); // site-free coupling
    arma::mat m(2, 2, arma::fill::eye);
    REQUIRE_NOTHROW(check_valid(Op("Matrix", int64_t(0), m)));
  }

  // --- Unknown type ---
  REQUIRE_THROWS_AS(check_valid(Op("NotAnOpType", int64_t(0))), Error);

  // --- site_required but no sites given ---
  REQUIRE_THROWS_AS(check_valid(Op("Sz")), Error);

  // --- wrong number of sites (Sz needs exactly 1) ---
  REQUIRE_THROWS_AS(
      check_valid(Op("Sz", std::vector<int64_t>{0, 1})), Error);

  // --- non-distinct sites where overlapping is disallowed
  //     (ScalarChirality needs 3 distinct sites) ---
  REQUIRE_THROWS_AS(
      check_valid(Op("ScalarChirality", std::vector<int64_t>{0, 0, 1})), Error);

  // --- matrix required but absent (Matrix type without a matrix) ---
  REQUIRE_THROWS_AS(check_valid(Op("Matrix", int64_t(0))), Error);

  // --- matrix present but not allowed (Sz carrying a matrix) ---
  {
    arma::mat m(2, 2, arma::fill::eye);
    REQUIRE_THROWS_AS(check_valid(Op("Sz", int64_t(0), m)), Error);
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
}

TEST_CASE("check_valid Monomial and OpSum", "[operators]") try {
  using namespace xdiag;
  using namespace xdiag::operators;

  Op good("Sz", int64_t(0));
  Op bad("Sz"); // site_required but missing

  // --- Monomial ---
  REQUIRE_NOTHROW(check_valid(Monomial{good, good}));
  REQUIRE_THROWS_AS(check_valid(Monomial{good, bad}), Error);

  // --- OpSum ---
  {
    OpSum ops;
    ops += good;
    REQUIRE_NOTHROW(check_valid(ops));
  }
  {
    OpSum ops;
    ops += bad;
    REQUIRE_THROWS_AS(check_valid(ops), Error);
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
}

TEST_CASE("valid must_* helpers", "[operators]") try {
  using namespace xdiag;
  using namespace xdiag::operators;

  Op sited("Sz", int64_t(2));
  Op sitefree("HubbardU");
  Op two("SdotS", std::vector<int64_t>{0, 1});
  Op dup("SdotS", std::vector<int64_t>{1, 1});
  arma::mat m(2, 2, arma::fill::eye);
  Op withmat("Matrix", int64_t(0), m);

  // --- must_have_sites ---
  REQUIRE_NOTHROW(must_have_sites(sited));
  REQUIRE_THROWS_AS(must_have_sites(sitefree), Error);

  // --- must_not_have_sites ---
  REQUIRE_NOTHROW(must_not_have_sites(sitefree));
  REQUIRE_THROWS_AS(must_not_have_sites(sited), Error);

  // --- must_have_nsites ---
  REQUIRE_NOTHROW(must_have_nsites(two, 2));
  REQUIRE_THROWS_AS(must_have_nsites(two, 3), Error);   // wrong count
  REQUIRE_THROWS_AS(must_have_nsites(sitefree, 1), Error); // no sites at all

  // --- must_have_disjoint_sites ---
  REQUIRE_NOTHROW(must_have_disjoint_sites(two));
  REQUIRE_THROWS_AS(must_have_disjoint_sites(dup), Error);      // duplicate
  REQUIRE_THROWS_AS(must_have_disjoint_sites(sitefree), Error); // no sites

  // --- must_have_sites_in_range: valid range [0, 4) ---
  REQUIRE_NOTHROW(must_have_sites_in_range(sited, 0, 4));
  REQUIRE_THROWS_AS(must_have_sites_in_range(sited, 0, 2), Error); // site 2 >= 2
  REQUIRE_THROWS_AS(must_have_sites_in_range(sited, 3, 5), Error); // site 2 < 3
  REQUIRE_THROWS_AS(must_have_sites_in_range(sitefree, 0, 4), Error); // no sites

  // --- must_have_matrix / must_not_have_matrix ---
  REQUIRE_NOTHROW(must_have_matrix(withmat));
  REQUIRE_THROWS_AS(must_have_matrix(sited), Error);
  REQUIRE_NOTHROW(must_not_have_matrix(sited));
  REQUIRE_THROWS_AS(must_not_have_matrix(withmat), Error);

} catch (xdiag::Error const &e) {
  error_trace(e);
}
