// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/algebra/utils/combine_matrix_ops.hpp>
#include <xdiag/utils/logger.hpp>

// Indexing convention reminder:
//   For an Op on sites [s0, s1, ...], the matrix row/column index k encodes
//   the local spin states as k = σ_{s0}*1 + σ_{s1}*2 + σ_{s2}*4 + ...
//   (bit j of k = spin at s_j).
//
//   In armadillo's kron(A, B): B is the inner (low-bit) factor, A the outer
//   (high-bit) factor. So for sites [s0, s1], a matrix M where s0 is inner
//   and s1 is outer is written arma::kron(M_s1, M_s0).
//
//   Consequence for combine_matrix_ops with two disjoint ops on disjoint
//   sites [s0_list, s1_list]:
//     combined_matrix = arma::kron(M1, M0)
//   (M0 is inner because s0_list appear first = lower bits)

TEST_CASE("combine_matrix_ops", "[operators]") try {
  using namespace xdiag;
  using namespace xdiag::algebra;
  using namespace arma;

  Log("Testing combine_matrix_ops");

  // Standard 1-site spin-1/2 matrices
  mat sz({{0.5, 0.0}, {0.0, -0.5}});
  mat sx({{0.0, 0.5}, {0.5, 0.0}});
  cx_mat sy(mat({{0., 0.}, {0., 0.}}), mat({{0., -0.5}, {0.5, 0.}}));
  mat id2 = eye(2, 2);

  // --- Single op: combining one op returns an equivalent op
  {
    mat M = kron(sx, sz); // sx on bit1=site5, sz on bit0=site2
    auto op = Op("Matrix", std::vector<int64_t>{2, 5}, M);
    auto result = combine_matrix_ops({op}, 2);
    REQUIRE(result.sites() == std::vector<int64_t>({2, 5}));
    REQUIRE(norm(result.matrix().as<arma::mat>() - M) < 1e-12);
  }

  // --- Two disjoint single-site ops, both real
  // op0 = Sz on site 0, op1 = Sx on site 1
  // all_sites = [0, 1]
  // embed(op0) = kron(id2, sz), embed(op1) = kron(sx, id2)
  // combined   = kron(id2,sz) * kron(sx,id2) = kron(sx, sz)
  {
    auto result =
        combine_matrix_ops({Op("Matrix", 0, sz), Op("Matrix", 1, sx)}, 2);
    REQUIRE(result.sites() == std::vector<int64_t>({0, 1}));
    REQUIRE(norm(result.matrix().as<arma::mat>() - kron(sx, sz)) < 1e-12);
  }

  // --- Two disjoint single-site ops, both complex
  // op0 = Sy on site 0, op1 = Sy on site 1
  // combined = kron(sy, sy)
  {
    auto result =
        combine_matrix_ops({Op("Matrix", 0, sy), Op("Matrix", 1, sy)}, 2);
    REQUIRE(result.sites() == std::vector<int64_t>({0, 1}));
    REQUIRE(norm(result.matrix().as<arma::cx_mat>() - kron(sy, sy)) < 1e-12);
  }

  // --- Mixed real/complex: one real op, one complex op
  // op0 = Sz on site 0 (real), op1 = Sy on site 1 (complex)
  // combined = kron(sy, cx_mat(sz))
  {
    auto result =
        combine_matrix_ops({Op("Matrix", 0, sz), Op("Matrix", 1, sy)}, 2);
    REQUIRE(result.sites() == std::vector<int64_t>({0, 1}));
    cx_mat expected = kron(sy, conv_to<cx_mat>::from(sz));
    REQUIRE(norm(result.matrix().as<arma::cx_mat>() - expected) < 1e-12);
  }

  // --- Two disjoint multi-site ops
  // op0 on sites [0, 1] with M0 = kron(sx, sz)
  // op1 on sites [2, 3] with M1 = kron(sz, sx)
  // all_sites = [0, 1, 2, 3]
  // embed(op0) = kron(id4, M0), embed(op1) = kron(M1, id4)
  // combined   = kron(M1, M0)
  {
    mat M0 = kron(sx, sz);
    mat M1 = kron(sz, sx);
    auto result =
        combine_matrix_ops({Op("Matrix", std::vector<int64_t>{0, 1}, M0),
                            Op("Matrix", std::vector<int64_t>{2, 3}, M1)}, 2);
    REQUIRE(result.sites() == std::vector<int64_t>({0, 1, 2, 3}));
    REQUIRE(norm(result.matrix().as<arma::mat>() - kron(M1, M0)) < 1e-12);
  }

  // --- Three disjoint single-site ops
  // op0 = Sz on site 0, op1 = Sx on site 1, op2 = Sz on site 2
  // combined = embed(op2) * embed(op1) * embed(op0) ... no, combined =
  //   embed(op0) * embed(op1) * embed(op2) where each embed acts on [0,1,2]:
  //   embed(op0) = kron(id2, kron(id2, sz))
  //   embed(op1) = kron(id2, kron(sx, id2))
  //   embed(op2) = kron(sz, kron(id2, id2))
  //   combined   = kron(sz, kron(sx, sz))
  {
    auto result = combine_matrix_ops({Op("Matrix", 0, sz), Op("Matrix", 1, sx),
                                      Op("Matrix", 2, sz)}, 2);
    REQUIRE(result.sites() == std::vector<int64_t>({0, 1, 2}));
    mat expected = kron(sz, kron(sx, sz));
    REQUIRE(norm(result.matrix().as<arma::mat>() - expected) < 1e-12);
  }

  // --- Overlapping sites: same single site, matrix product
  // op0 = Sz on site 0, op1 = Sx on site 0
  // all_sites = [0], combined = Sz * Sx (ordinary matrix multiplication)
  {
    auto result =
        combine_matrix_ops({Op("Matrix", 0, sz), Op("Matrix", 0, sx)}, 2);
    REQUIRE(result.sites() == std::vector<int64_t>({0}));
    REQUIRE(norm(result.matrix().as<arma::mat>() - sz * sx) < 1e-12);
  }

  // --- Overlapping sites: same site composed three times
  // Sz * Sx * Sz on site 0
  {
    auto result = combine_matrix_ops(
        {Op("Matrix", 0, sz), Op("Matrix", 0, sx), Op("Matrix", 0, sz)}, 2);
    REQUIRE(result.sites() == std::vector<int64_t>({0}));
    REQUIRE(norm(result.matrix().as<arma::mat>() - sz * sx * sz) < 1e-12);
  }

  // --- Overlapping sites: partial overlap between two 1-site ops and a
  //     second 1-site op on a fresh site
  // op0 = Sz on site 0, op1 = Sx on site 0, op2 = Sy on site 1
  // all_sites = [0, 1]
  // embed(op0) = kron(id2, sz), embed(op1) = kron(id2, sx),
  // embed(op2) = kron(sy, id2)
  // combined   = kron(id2,sz)*kron(id2,sx)*kron(sy,id2)
  //            = kron(id2, sz*sx) * kron(sy, id2)
  //            = kron(sy, sz*sx)
  {
    auto result = combine_matrix_ops({Op("Matrix", 0, sz), Op("Matrix", 0, sx),
                                      Op("Matrix", 1, sy)}, 2);
    REQUIRE(result.sites() == std::vector<int64_t>({0, 1}));
    cx_mat expected = kron(sy, conv_to<cx_mat>::from(sz * sx));
    REQUIRE(norm(result.matrix().as<arma::cx_mat>() - expected) < 1e-12);
  }

  // --- Overlapping sites: two 2-site ops sharing one site
  // op0 on sites [0, 1] with M0 = kron(sz, sz)  (Sz⊗Sz)
  // op1 on site  [1]    with M1 = sz             (Sz)
  // all_sites = [0, 1]
  // embed(op0) = kron(sz, sz)      (covers all sites exactly)
  // embed(op1) = kron(sz, id2)     (Sz on bit1=site1, I on bit0=site0)
  // combined   = kron(sz,sz)*kron(sz,id2) = kron(sz*sz, sz*id2)
  //            = kron(0.25*id2, sz) = 0.25 * kron(id2, sz)
  {
    mat M0 = kron(sz, sz);
    auto result =
        combine_matrix_ops({Op("Matrix", std::vector<int64_t>{0, 1}, M0),
                            Op("Matrix", 1, sz)}, 2);
    REQUIRE(result.sites() == std::vector<int64_t>({0, 1}));
    mat expected = 0.25 * kron(id2, sz);
    REQUIRE(norm(result.matrix().as<arma::mat>() - expected) < 1e-12);
  }

  // --- Sites in non-contiguous, non-ascending order
  // op0 = Sz on site 5, op1 = Sx on site 2
  // all_sites = [5, 2] (first appearance order)
  // embed(op0) = kron(id2, sz) (sz on bit0=site5 inner, I on bit1=site2 outer)
  // embed(op1) = kron(sx, id2) (I on bit0=site5, sx on bit1=site2 outer)
  // combined   = kron(sx, sz)
  {
    auto result =
        combine_matrix_ops({Op("Matrix", 5, sz), Op("Matrix", 2, sx)}, 2);
    REQUIRE(result.sites() == std::vector<int64_t>({5, 2}));
    REQUIRE(norm(result.matrix().as<arma::mat>() - kron(sx, sz)) < 1e-12);
  }

  Log("Done testing combine_matrix_ops");
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
  throw;
}
