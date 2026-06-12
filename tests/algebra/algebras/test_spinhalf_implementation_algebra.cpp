// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/algebra/algebras/matrix_algebra.hpp>
#include <xdiag/algebra/algebras/spinhalf_implementation_algebra.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace xdiag::algebra;
using namespace arma;

TEST_CASE("spinhalf_implementation_algebra", "[operators]") try {

  auto impl_alg = spinhalf_implementation_algebra(4);
  auto mat_alg = matrix_algebra(4, 2);

  // Standard spin-1/2 matrices (for direct matrix comparisons). Local basis
  // index 0 = down (m = -1/2), index 1 = up (m = +1/2), matching the basis bit
  // convention (a set bit is "up"); see op_to_matrix_op.
  mat sp = {{0.0, 0.0}, {1.0, 0.0}};
  mat sm = {{0.0, 1.0}, {0.0, 0.0}};
  mat sz = {{-0.5, 0.0}, {0.0, 0.5}};
  mat sx = {{0.0, 0.5}, {0.5, 0.0}};
  cx_mat sy(mat({{0., 0.}, {0., 0.}}), mat({{0., 0.5}, {-0.5, 0.}}));
  mat id2 = eye(2, 2);

  // -------------------------------------------------------------------------
  // Protected single-site operators: Sz, S+, S- stay as named ops
  // -------------------------------------------------------------------------
  Log("Testing spinhalf_implementation_algebra: protected single-site "
      "operators");
  {
    auto r = normal_order(OpSum(Op("Sz", 0)), impl_alg);
    REQUIRE(r.size() == 1);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Sz");
      REQUIRE(mono[0].sites() == std::vector<int64_t>{0});
    }
    REQUIRE(isapprox(r, OpSum(Op("Sz", 0)), mat_alg));
  }
  {
    auto r = normal_order(OpSum(Op("S+", 0)), impl_alg);
    REQUIRE(r.size() == 1);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "S+");
    }
    REQUIRE(isapprox(r, OpSum(Op("S+", 0)), mat_alg));
  }
  {
    auto r = normal_order(OpSum(Op("S-", 0)), impl_alg);
    REQUIRE(r.size() == 1);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "S-");
    }
    REQUIRE(isapprox(r, OpSum(Op("S-", 0)), mat_alg));
  }

  // -------------------------------------------------------------------------
  // Protected two-site operators: SzSz, Exchange, ScalarChirality stay named
  // -------------------------------------------------------------------------
  Log("Testing spinhalf_implementation_algebra: protected two-site operators");
  {
    auto r = normal_order(OpSum(Op("SzSz", {0, 1})), impl_alg);
    REQUIRE(r.size() == 1);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "SzSz");
      REQUIRE(mono[0].sites() == std::vector<int64_t>{0, 1});
    }
    REQUIRE(isapprox(r, OpSum(Op("SzSz", {0, 1})), mat_alg));
  }
  {
    auto r = normal_order(OpSum(Op("Exchange", {0, 1})), impl_alg);
    REQUIRE(r.size() == 1);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Exchange");
      REQUIRE(mono[0].sites() == std::vector<int64_t>{0, 1});
    }
    REQUIRE(isapprox(r, OpSum(Op("Exchange", {0, 1})), mat_alg));
  }
  {
    auto r = normal_order(OpSum(Op("ScalarChirality", {0, 1, 2})), impl_alg);
    REQUIRE(r.size() == 1);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "ScalarChirality");
      REQUIRE(mono[0].sites() == std::vector<int64_t>{0, 1, 2});
    }
  }

  // -------------------------------------------------------------------------
  // Id and a pre-built Matrix op stay as-is
  // -------------------------------------------------------------------------
  Log("Testing spinhalf_implementation_algebra: Id and Matrix pass-through");
  {
    auto r = normal_order(OpSum(Op("Id")), impl_alg);
    REQUIRE(r.size() == 1);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Id");
    }
    REQUIRE(isapprox(r, OpSum(Op("Id")), mat_alg));
  }
  {
    auto r = normal_order(OpSum(Op("Matrix", 0, sz)), impl_alg);
    REQUIRE(r.size() == 1);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Matrix");
    }
    REQUIRE(isapprox(r, OpSum(Op("Matrix", 0, sz)), mat_alg));
  }

  // -------------------------------------------------------------------------
  // Sx and Sy are NOT protected: they become Matrix ops
  // -------------------------------------------------------------------------
  Log("Testing spinhalf_implementation_algebra: Sx and Sy become Matrix");
  {
    auto r = normal_order(OpSum(Op("Sx", 0)), impl_alg);
    REQUIRE(r.size() == 1);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Matrix");
    }
    REQUIRE(isapprox(r, OpSum(Op("Matrix", 0, sx)), mat_alg));
  }
  {
    auto r = normal_order(OpSum(Op("Sy", 0)), impl_alg);
    REQUIRE(r.size() == 1);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Matrix");
    }
    REQUIRE(isapprox(r, OpSum(Op("Matrix", 0, sy)), mat_alg));
  }
  {
    // Sx{0} and Sy{0} give the same result as in matrix_algebra
    REQUIRE(isapprox(normal_order(OpSum(Op("Sx", 0)), impl_alg),
                     normal_order(OpSum(Op("Sx", 0)), mat_alg), mat_alg));
    REQUIRE(isapprox(normal_order(OpSum(Op("Sy", 0)), impl_alg),
                     normal_order(OpSum(Op("Sy", 0)), mat_alg), mat_alg));
  }

  // -------------------------------------------------------------------------
  // SdotS{i,j} → Exchange{i,j} + SzSz{i,j}  (not to S+/S-/Sz or Matrix)
  // -------------------------------------------------------------------------
  Log("Testing spinhalf_implementation_algebra: SdotS expansion");
  {
    auto r = normal_order(OpSum(Op("SdotS", {0, 1})), impl_alg);
    REQUIRE(r.size() == 2);
    bool has_exchange = false, has_szsz = false;
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      if (mono[0].type() == "Exchange") {
        has_exchange = true;
        REQUIRE(mono[0].sites() == std::vector<int64_t>{0, 1});
      } else if (mono[0].type() == "SzSz") {
        has_szsz = true;
        REQUIRE(mono[0].sites() == std::vector<int64_t>{0, 1});
      } else {
        FAIL("SdotS must expand only to Exchange and SzSz, got: " +
             mono[0].type());
      }
    }
    REQUIRE(has_exchange);
    REQUIRE(has_szsz);
    // Equivalent to explicit sum
    REQUIRE(
        isapprox(r, OpSum(Op("Exchange", {0, 1})) + OpSum(Op("SzSz", {0, 1})), mat_alg));
  }
  {
    // SdotS{0,0} same-site: still expands to Exchange{0,0} + SzSz{0,0}
    // (the named ops remain; semantically = 0.75 * Id for spin-1/2)
    auto r = normal_order(OpSum(Op("SdotS", {0, 0})), impl_alg);
    REQUIRE(r.size() == 2);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      auto t = mono[0].type();
      REQUIRE((t == "Exchange" || t == "SzSz"));
    }
  }
  {
    // Coefficient is preserved: 2.0 * SdotS{0,1} → 2*Exchange + 2*SzSz
    auto r =
        normal_order(2.0 * Op("SdotS", std::vector<int64_t>{0, 1}), impl_alg);
    for (auto const &[c, mono] : r) {
      REQUIRE(std::abs(c.scalar().real() - 2.0) < 1e-12);
    }
  }

  // -------------------------------------------------------------------------
  // Multi-op monomials collapse to a single Matrix op and match matrix_algebra
  // -------------------------------------------------------------------------
  Log("Testing spinhalf_implementation_algebra: multi-op monomials → Matrix");
  {
    // Sz{0}*Sz{1}
    auto r_impl = normal_order(1.0 * (Op("Sz", 0) * Op("Sz", 1)), impl_alg);
    auto r_mat = normal_order(1.0 * (Op("Sz", 0) * Op("Sz", 1)), mat_alg);
    for (auto const &[c, mono] : r_impl) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Matrix");
    }
    REQUIRE(isapprox(r_impl, r_mat, mat_alg));
  }
  {
    // S+{0}*S-{1}
    auto r_impl = normal_order(1.0 * (Op("S+", 0) * Op("S-", 1)), impl_alg);
    auto r_mat = normal_order(1.0 * (Op("S+", 0) * Op("S-", 1)), mat_alg);
    REQUIRE(isapprox(r_impl, r_mat, mat_alg));
  }
  {
    // S-{0}*S+{1}
    auto r_impl = normal_order(1.0 * (Op("S-", 0) * Op("S+", 1)), impl_alg);
    auto r_mat = normal_order(1.0 * (Op("S-", 0) * Op("S+", 1)), mat_alg);
    REQUIRE(isapprox(r_impl, r_mat, mat_alg));
  }
  {
    // S+{0}*S+{0} same-site → zero matrix
    auto r_impl = normal_order(1.0 * (Op("S+", 0) * Op("S+", 0)), impl_alg);
    auto r_mat = normal_order(1.0 * (Op("S+", 0) * Op("S+", 0)), mat_alg);
    REQUIRE(isapprox(r_impl, r_mat, mat_alg));
  }
  {
    // Sx{0}*Sz{1}: Sx is not protected → first becomes Matrix, then combined
    auto r_impl = normal_order(1.0 * (Op("Sx", 0) * Op("Sz", 1)), impl_alg);
    auto r_mat = normal_order(1.0 * (Op("Sx", 0) * Op("Sz", 1)), mat_alg);
    REQUIRE(isapprox(r_impl, r_mat, mat_alg));
  }
  {
    // Three-site: Sz{0}*S+{1}*S-{2}
    auto r_impl =
        normal_order(1.0 * (Op("Sz", 0) * Op("S+", 1) * Op("S-", 2)), impl_alg);
    auto r_mat =
        normal_order(1.0 * (Op("Sz", 0) * Op("S+", 1) * Op("S-", 2)), mat_alg);
    for (auto const &[c, mono] : r_impl) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Matrix");
      REQUIRE(mono[0].sites() == std::vector<int64_t>({0, 1, 2}));
    }
    REQUIRE(isapprox(r_impl, r_mat, mat_alg));
  }
  {
    // Non-adjacent sites: Sz{3}*S+{1} — result has sorted sites
    auto r = normal_order(1.0 * (Op("Sz", 3) * Op("S+", 1)), impl_alg);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Matrix");
      auto const &sites = mono[0].sites();
      for (int64_t i = 0; i + 1 < (int64_t)sites.size(); ++i)
        REQUIRE(sites[i] < sites[i + 1]);
    }
    REQUIRE(
        isapprox(r, normal_order(1.0 * (Op("Sz", 3) * Op("S+", 1)), mat_alg), mat_alg));
  }

  // -------------------------------------------------------------------------
  // Heisenberg chain: SdotS terms expand to Exchange + SzSz
  // -------------------------------------------------------------------------
  Log("Testing spinhalf_implementation_algebra: Heisenberg chain");
  {
    // H = SdotS{0,1} + SdotS{1,2}  →  Exchange{0,1} + SzSz{0,1}
    //                                 + Exchange{1,2} + SzSz{1,2}
    auto H = OpSum(Op("SdotS", {0, 1})) + OpSum(Op("SdotS", {1, 2}));
    auto r = normal_order(H, impl_alg);
    REQUIRE(r.size() == 4);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      auto t = mono[0].type();
      REQUIRE((t == "Exchange" || t == "SzSz"));
    }
    REQUIRE(isapprox(
        r, OpSum(Op("Exchange", {0, 1})) + OpSum(Op("SzSz", {0, 1})) +
               OpSum(Op("Exchange", {1, 2})) + OpSum(Op("SzSz", {1, 2})), mat_alg));
  }
  {
    // Mix of protected and non-protected: J*SdotS + h*Sz
    auto H =
        OpSum(Op("SdotS", {0, 1})) + OpSum(Op("Sz", 0)) + OpSum(Op("Sz", 1));
    auto r = normal_order(H, impl_alg);
    bool has_exchange = false, has_szsz = false, has_sz = false;
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      if (mono[0].type() == "Exchange")
        has_exchange = true;
      else if (mono[0].type() == "SzSz")
        has_szsz = true;
      else if (mono[0].type() == "Sz")
        has_sz = true;
      else
        FAIL("unexpected type: " + mono[0].type());
    }
    REQUIRE(has_exchange);
    REQUIRE(has_szsz);
    REQUIRE(has_sz);
  }

  // -------------------------------------------------------------------------
  // Consistency with matrix_algebra: Exchange and SzSz from spinhalf_impl
  // produce the same physical matrix when evaluated via matrix_algebra
  // -------------------------------------------------------------------------
  Log("Testing spinhalf_implementation_algebra: consistency with "
      "matrix_algebra");
  {
    // SdotS{0,1} via matrix_algebra gives ONE Matrix op (full Heisenberg 4×4).
    // SdotS{0,1} via spinhalf_impl gives Exchange + SzSz; running those
    // through matrix_algebra also gives ONE Matrix op (after collect merges
    // the two Matrix ops with the same sites into one).  The results must
    // agree.
    auto sdots_via_mat = normal_order(OpSum(Op("SdotS", {0, 1})), mat_alg);
    auto impl_result = normal_order(OpSum(Op("SdotS", {0, 1})), impl_alg);
    auto impl_via_mat = normal_order(impl_result, mat_alg);

    REQUIRE(sdots_via_mat.size() == 1);
    REQUIRE(impl_via_mat.size() == 1);
    REQUIRE(isapprox(sdots_via_mat, impl_via_mat, mat_alg));
  }

  // -------------------------------------------------------------------------
  // Allowed types validation
  // -------------------------------------------------------------------------
  Log("Testing spinhalf_implementation_algebra: allowed types validation");
  {
    // REQUIRE_THROWS(normal_order(OpSum(Op("Cdagup", 0)), impl_alg));
    // REQUIRE_THROWS(normal_order(OpSum(Op("Nup", 0)), impl_alg));
    // REQUIRE_THROWS(normal_order(OpSum(Op("Hop", {0, 1})), impl_alg));
  }

  Log("Done testing spinhalf_implementation_algebra");
} catch (xdiag::Error const &e) {
  error_trace(e);
}
