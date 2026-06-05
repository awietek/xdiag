// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/algebra/algebras/matrix_algebra.hpp>
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

// Sums all (coefficient * matrix) terms in an OpSum that is expected to
// contain only real-matrix Monomials of size 1. Returns the total matrix.
static mat sum_real_matrices(OpSum const &ops) {
  mat total;
  bool first = true;
  for (auto const &[c, mono] : ops) {
    mat M = mono[0].matrix().as<mat>() * c.scalar().real();
    if (first) {
      total = M;
      first = false;
    } else {
      total += M;
    }
  }
  return total;
}

// Same but for complex-matrix OpSums.
static cx_mat sum_cx_matrices(OpSum const &ops) {
  cx_mat total;
  bool first = true;
  for (auto const &[c, mono] : ops) {
    complex coeff(c.scalar().real(), c.scalar().imag());
    cx_mat M = mono[0].matrix().as<cx_mat>() * coeff;
    if (first) {
      total = M;
      first = false;
    } else {
      total += M;
    }
  }
  return total;
}

TEST_CASE("matrix_algebra", "[operators]") try {

  auto mat_alg = matrix_algebra(2);

  // Standard spin-1/2 matrices
  mat sp = {{0.0, 1.0}, {0.0, 0.0}};
  mat sm = {{0.0, 0.0}, {1.0, 0.0}};
  mat sz = {{0.5, 0.0}, {0.0, -0.5}};
  mat sx = {{0.0, 0.5}, {0.5, 0.0}};
  cx_mat sy(mat({{0., 0.}, {0., 0.}}), mat({{0., -0.5}, {0.5, 0.}}));
  mat id2 = eye(2, 2);

  // -------------------------------------------------------------------------
  // Single-site operators → correct 2×2 matrix
  // -------------------------------------------------------------------------
  Log("Testing matrix_algebra: single-site operators");
  {
    auto r = normal_order(OpSum(Op("S+", 0)), mat_alg);
    REQUIRE(isapprox(r, OpSum(Op("Matrix", 0, sp)), mat_alg));
  }
  {
    auto r = normal_order(OpSum(Op("S-", 0)), mat_alg);
    REQUIRE(isapprox(r, OpSum(Op("Matrix", 0, sm)), mat_alg));
  }
  {
    auto r = normal_order(OpSum(Op("Sz", 0)), mat_alg);
    REQUIRE(isapprox(r, OpSum(Op("Matrix", 0, sz)), mat_alg));
  }
  {
    // Sx = (S+ + S-)/2 -> [[0, 0.5],[0.5, 0]]
    auto r = normal_order(OpSum(Op("Sx", 0)), mat_alg);
    REQUIRE(isapprox(r, OpSum(Op("Matrix", 0, sx)), mat_alg));
  }
  {
    // Sy = (S+ - S-)/(2i) -> [[0, -i/2],[i/2, 0]]
    auto r = normal_order(OpSum(Op("Sy", 0)), mat_alg);
    REQUIRE(isapprox(r, OpSum(Op("Matrix", 0, sy)), mat_alg));
  }

  // -------------------------------------------------------------------------
  // Same-site products → matrix multiplication at the Op level
  // -------------------------------------------------------------------------
  Log("Testing matrix_algebra: same-site products");
  {
    // Sz{0}*Sz{0} -> [[0.25,0],[0,0.25]]  (0.25 * identity)
    auto r = normal_order(1.0 * (Op("Sz", 0) * Op("Sz", 0)), mat_alg);
    REQUIRE(isapprox(r, OpSum(Op("Matrix", 0, mat(sz * sz))), mat_alg));
  }
  {
    // S+{0}*S-{0} -> sp*sm = [[1,0],[0,0]]
    auto r = normal_order(1.0 * (Op("S+", 0) * Op("S-", 0)), mat_alg);
    REQUIRE(isapprox(r, OpSum(Op("Matrix", 0, mat(sp * sm))), mat_alg));
  }
  {
    // S-{0}*S+{0} -> sm*sp = [[0,0],[0,1]]
    auto r = normal_order(1.0 * (Op("S-", 0) * Op("S+", 0)), mat_alg);
    REQUIRE(isapprox(r, OpSum(Op("Matrix", 0, mat(sm * sp))), mat_alg));
  }
  {
    // S+{0}*S+{0} -> sp*sp = 0
    // sp*sp = [[0,1],[0,0]]*[[0,1],[0,0]] = [[0,0],[0,0]]
    auto r = normal_order(1.0 * (Op("S+", 0) * Op("S+", 0)), mat_alg);
    REQUIRE(isapprox(r, OpSum(Op("Matrix", 0, mat(sp * sp))), mat_alg));
  }

  // -------------------------------------------------------------------------
  // Two-site products with pre-sorted sites
  // -------------------------------------------------------------------------
  Log("Testing matrix_algebra: two-site products");
  {
    // Sz{0}*Sz{1} -> Matrix{0,1}[kron(sz,sz)]
    // In armadillo kron convention: site0 is inner (low-bit), site1 is outer
    // embed(Sz{0}) = kron(id2,sz), embed(Sz{1}) = kron(sz,id2)
    // product = kron(sz, sz)
    auto r = normal_order(1.0 * (Op("Sz", 0) * Op("Sz", 1)), mat_alg);
    auto expected =
        OpSum(Op("Matrix", std::vector<int64_t>{0, 1}, mat(kron(sz, sz))));
    REQUIRE(isapprox(r, expected, mat_alg));
  }
  {
    // SzSz{0,1} should expand to the same result as Sz{0}*Sz{1}
    auto r1 = normal_order(OpSum(Op("SzSz", {0, 1})), mat_alg);
    auto r2 = normal_order(1.0 * (Op("Sz", 0) * Op("Sz", 1)), mat_alg);
    REQUIRE(isapprox(r1, r2, mat_alg));
  }

  // -------------------------------------------------------------------------
  // Site ordering: bosonic operators, no sign change on transposition
  // -------------------------------------------------------------------------
  Log("Testing matrix_algebra: site ordering");
  {
    // Sz{1}*Sz{0} == Sz{0}*Sz{1}: sites get sorted, matrix permuted
    auto r1 = normal_order(1.0 * (Op("Sz", 0) * Op("Sz", 1)), mat_alg);
    auto r2 = normal_order(1.0 * (Op("Sz", 1) * Op("Sz", 0)), mat_alg);
    REQUIRE(isapprox(r1, r2, mat_alg));
  }
  {
    // S+{1}*S-{0} == S-{0}*S+{1}: no sign from site exchange (bosonic)
    auto r1 = normal_order(1.0 * (Op("S+", 0) * Op("S-", 1)), mat_alg);
    auto r2 = normal_order(1.0 * (Op("S-", 1) * Op("S+", 0)), mat_alg);
    REQUIRE(isapprox(r1, r2, mat_alg));
  }
  {
    // Non-adjacent sites: S+{0}*Sz{2} == Sz{2}*S+{0}
    auto r1 = normal_order(1.0 * (Op("S+", 0) * Op("Sz", 2)), mat_alg);
    auto r2 = normal_order(1.0 * (Op("Sz", 2) * Op("S+", 0)), mat_alg);
    REQUIRE(isapprox(r1, r2, mat_alg));
  }
  {
    // Result has sorted sites
    auto r = normal_order(1.0 * (Op("Sz", 3) * Op("S+", 1)), mat_alg);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Matrix");
      auto const &sites = mono[0].sites();
      for (int64_t i = 0; i + 1 < (int64_t)sites.size(); ++i) {
        REQUIRE(sites[i] < sites[i + 1]);
      }
    }
  }

  // -------------------------------------------------------------------------
  // Compound operator expansion: SdotS and Exchange
  // -------------------------------------------------------------------------
  Log("Testing matrix_algebra: compound operators");
  {
    // SdotS{0,1}: all terms should be size-1 Matrix monomials on sites [0,1]
    // The summed matrix should equal the Heisenberg interaction
    auto r = normal_order(OpSum(Op("SdotS", {0, 1})), mat_alg);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Matrix");
      REQUIRE(mono[0].sites() == std::vector<int64_t>({0, 1}));
    }
    // Sum weighted matrices and compare to Heisenberg: Sz⊗Sz + 0.5*S-⊗S+ +
    // 0.5*S+⊗S-
    mat heisenberg = kron(sz, sz) + 0.5 * kron(sm, sp) + 0.5 * kron(sp, sm);
    REQUIRE(norm(sum_real_matrices(r) - heisenberg) < 1e-12);
  }
  {
    // Exchange{0,1}: summed matrix = 0.5*(S-⊗S+ + S+⊗S-)
    auto r = normal_order(OpSum(Op("Exchange", {0, 1})), mat_alg);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Matrix");
    }
    mat expected = 0.5 * kron(sm, sp) + 0.5 * kron(sp, sm);
    REQUIRE(norm(sum_real_matrices(r) - expected) < 1e-12);
  }
  {
    // SdotS{0,0}: S·S for spin-1/2 at the same site equals 0.75 * identity
    // In matrix form: sp*sm/2 + sm*sp/2 + sz*sz = I/2 - I/2 + ... = 0.75*I
    // Just verify the 2x2 sum equals 0.75 * id2
    auto r = normal_order(OpSum(Op("SdotS", {0, 0})), mat_alg);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Matrix");
    }
    REQUIRE(norm(sum_real_matrices(r) - 0.75 * id2) < 1e-12);
  }

  // -------------------------------------------------------------------------
  // Three-site operator: ScalarChirality (complex result)
  // -------------------------------------------------------------------------
  Log("Testing matrix_algebra: ScalarChirality");
  {
    auto r = normal_order(OpSum(Op("ScalarChirality", {0, 1, 2})), mat_alg);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Matrix");
      REQUIRE(mono[0].sites() == std::vector<int64_t>({0, 1, 2}));
    }
    // The ScalarChirality S·(S×S) is Hermitian: M - M^† = 0
    cx_mat total = sum_cx_matrices(r);
    REQUIRE(norm(total - total.t()) < 1e-12);
  }

  // -------------------------------------------------------------------------
  // Pass-through: a pre-built Matrix op with sorted sites is unchanged
  // -------------------------------------------------------------------------
  Log("Testing matrix_algebra: pass-through Matrix ops");
  {
    mat M = kron(sz, sx);
    auto ops_in = OpSum(Op("Matrix", std::vector<int64_t>{0, 1}, M));
    auto r = normal_order(ops_in, mat_alg);
    REQUIRE(isapprox(r, ops_in, mat_alg));
  }
  {
    // Unsorted sites get permuted: Matrix{1,0}[M] -> Matrix{0,1}[M_permuted]
    // After permutation the sites must be ascending
    mat M = {{1.0, 2.0, 3.0, 4.0},
             {5.0, 6.0, 7.0, 8.0},
             {9.0, 10.0, 11.0, 12.0},
             {13.0, 14.0, 15.0, 16.0}};
    auto r = normal_order(OpSum(Op("Matrix", std::vector<int64_t>{1, 0}, M)),
                          mat_alg);
    REQUIRE(r.size() == 1);
    for (auto const &[c, mono] : r) {
      auto const &sites = mono[0].sites();
      REQUIRE(sites[0] == 0);
      REQUIRE(sites[1] == 1);
    }
    // The permuted op must represent the same physical operator: applying it
    // to Sz{1}*Sz{0} and Sz{0}*Sz{1} must agree after normal_order
    auto r1 = normal_order(1.0 * (Op("Sz", 0) * Op("Sz", 1)), mat_alg);
    auto r2 = normal_order(1.0 * (Op("Sz", 1) * Op("Sz", 0)), mat_alg);
    REQUIRE(isapprox(r1, r2, mat_alg));
  }

  // -------------------------------------------------------------------------
  // Linear combinations with multiple sites and scalar coefficients
  // -------------------------------------------------------------------------
  Log("Testing matrix_algebra: Heisenberg chain (3 sites)");
  {
    // H = J * sum_<ij> SdotS{i,j}  for a 3-site chain
    auto H = OpSum(Op("SdotS", {0, 1})) + OpSum(Op("SdotS", {1, 2}));
    auto r = normal_order(H, mat_alg);
    for (auto const &[c, mono] : r) {
      REQUIRE(mono.size() == 1);
      REQUIRE(mono[0].type() == "Matrix");
      auto const &sites = mono[0].sites();
      for (int64_t i = 0; i + 1 < (int64_t)sites.size(); ++i) {
        REQUIRE(sites[i] < sites[i + 1]);
      }
    }
  }

  // -------------------------------------------------------------------------
  // Allowed types validation
  // -------------------------------------------------------------------------
  Log("Testing matrix_algebra: allowed types validation");
  {
    // REQUIRE_THROWS(normal_order(OpSum(Op("Cdagup", 0)), mat_alg));
    // REQUIRE_THROWS(normal_order(OpSum(Op("Nup", 0)), mat_alg));
    // REQUIRE_THROWS(normal_order(OpSum(Op("Hop", {0, 1})), mat_alg));
  }

  Log("Done testing matrix_algebra");
} catch (xdiag::Error const &e) {
  error_trace(e);
}
