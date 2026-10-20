// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <vector>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/operators/hc.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

// Checks conj, transpose and hc of each OpSum against its dense matrix in block.
template <class block_t>
void check_transpose_conj_hc_matrices(block_t const &block,
                                      std::vector<xdiag::OpSum> const &ops_list) {
  using namespace xdiag;
  for (auto const &ops : ops_list) {
    auto M = matrixC(ops, block);

    // conj: complex-conjugate the coefficient and every matrix entry
    REQUIRE(arma::approx_equal(matrixC(conj(ops), block), arma::conj(M),
                               "absdiff", 1e-12));
    // transpose: transpose, without conjugation
    REQUIRE(arma::approx_equal(matrixC(transpose(ops), block), M.st(), "absdiff",
                               1e-12));
    // hc: conjugate transpose
    REQUIRE(arma::approx_equal(matrixC(hc(ops), block), M.t(), "absdiff",
                               1e-12));
    // hc = conj(transpose(.)) = transpose(conj(.))
    REQUIRE(arma::approx_equal(matrixC(hc(ops), block),
                               matrixC(conj(transpose(ops)), block), "absdiff",
                               1e-12));
    REQUIRE(arma::approx_equal(matrixC(hc(ops), block),
                               matrixC(transpose(conj(ops)), block), "absdiff",
                               1e-12));
  }
}

TEST_CASE("transpose, conj and hc", "[operators]") try {
  using namespace xdiag;

  Log("Testing transpose, conj and hc");
  // --- Named ops, compared with their dense matrices on 3 spins ---
  // The complex coefficient makes conj act on the coefficient as well as on
  // the matrix.
  complex c(0.7, -0.4);
  check_transpose_conj_hc_matrices(
      Spinhalf(3),
      {
          OpSum(c, Op("Sx", 0)),
          OpSum(c, Op("Sy", 1)),
          OpSum(c, Op("Sz", 2)),
          OpSum(c, Op("S+", 0)),
          OpSum(c, Op("S-", 2)),
          OpSum(c, Op("SzSz", {0, 1})),
          OpSum(c, Op("SdotS", {0, 1})),
          OpSum(c, Op("Exchange", {0, 2})),
          OpSum(c, Op("ExchangeAsym", {1, 2})),
          OpSum(c, Op("ScalarChirality", {0, 1, 2})),
      });

  // --- Signs and partners for single named ops ---
  // Sy and ScalarChirality are imaginary: conj flips the sign.
  REQUIRE(conj(OpSum(Op("Sy", 0))) == OpSum(-1.0, Op("Sy", 0)));
  REQUIRE(transpose(OpSum(Op("Sy", 0))) == OpSum(-1.0, Op("Sy", 0)));
  REQUIRE(hc(OpSum(Op("Sy", 0))) == OpSum(Op("Sy", 0)));

  // ExchangeAsym is real and antisymmetric: transpose flips the sign, conj keeps
  // it.
  REQUIRE(conj(OpSum(Op("ExchangeAsym", {0, 1}))) ==
          OpSum(Op("ExchangeAsym", {0, 1})));
  REQUIRE(transpose(OpSum(Op("ExchangeAsym", {0, 1}))) ==
          OpSum(-1.0, Op("ExchangeAsym", {0, 1})));
  REQUIRE(hc(OpSum(Op("ExchangeAsym", {0, 1}))) ==
          OpSum(-1.0, Op("ExchangeAsym", {0, 1})));

  // Sx is real and symmetric: both leave it unchanged.
  REQUIRE(conj(OpSum(Op("Sx", 0))) == OpSum(Op("Sx", 0)));
  REQUIRE(transpose(OpSum(Op("Sx", 0))) == OpSum(Op("Sx", 0)));

  // S+ transposes to S-, and conj leaves it unchanged.
  REQUIRE(transpose(OpSum(Op("S+", 0))) == OpSum(Op("S-", 0)));
  REQUIRE(conj(OpSum(Op("S+", 0))) == OpSum(Op("S+", 0)));
  REQUIRE(hc(OpSum(Op("S+", 0))) == OpSum(Op("S-", 0)));

  // --- Monomials: transpose reverses the product, conj keeps the order ---
  Monomial m{Op("S+", 0), Op("S-", 1)};
  REQUIRE(transpose(m) == OpSum(Monomial{Op("S+", 1), Op("S-", 0)}));
  REQUIRE(conj(m) == OpSum(m));
  REQUIRE(hc(m) == OpSum(Monomial{Op("S+", 1), Op("S-", 0)}));

  // --- Matrix ops: the stored matrix is transposed or conjugated ---
  auto block1 = Spinhalf(1);
  arma::cx_mat A = arma::randn<arma::cx_mat>(2, 2);
  OpSum mat_ops(Op("Matrix", 0, A));
  REQUIRE(arma::approx_equal(matrixC(transpose(mat_ops), block1), A.st(),
                             "absdiff", 1e-12));
  REQUIRE(arma::approx_equal(matrixC(conj(mat_ops), block1), arma::conj(A),
                             "absdiff", 1e-12));
  REQUIRE(arma::approx_equal(matrixC(hc(mat_ops), block1), A.t(), "absdiff",
                             1e-12));

  // --- adjoint is the same map as hc ---
  REQUIRE(adjoint(OpSum(Op("S+", 0))) == hc(OpSum(Op("S+", 0))));
  REQUIRE(adjoint(OpSum(Op("ExchangeAsym", {0, 1}))) ==
          hc(OpSum(Op("ExchangeAsym", {0, 1}))));
  REQUIRE(adjoint(m) == hc(m));
  REQUIRE(adjoint(Op("S+", 0)) == hc(Op("S+", 0)));

  // --- Unknown op types are rejected ---
  REQUIRE_THROWS(transpose(OpSum(Op("NoSuchOp", 0))));
  REQUIRE_THROWS(conj(OpSum(Op("NoSuchOp", 0))));
  REQUIRE_THROWS(hc(OpSum(Op("NoSuchOp", 0))));

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("transpose, conj and hc for fermion, electron and tJ blocks",
          "[operators]") try {
  using namespace xdiag;

  // Full Hilbert spaces (no sector restriction), so Cdag and C both act within
  // the dense matrices. The spin-down signs of electrons and tJ depend on the
  // spin-up occupation, which the transpose must reproduce.
  complex c(0.6, 0.8);

  // Spinless fermions. Hop and HopAsym are dispatched directly.
  check_transpose_conj_hc_matrices(
      Fermion(3),
      {
          OpSum(c, Op("Cdag", 1)),
          OpSum(c, Op("C", 0)),
          OpSum(c, Op("Hop", {0, 2})),
          OpSum(c, Op("HopAsym", {0, 2})),
          OpSum(c, Monomial{Op("Cdag", 0), Op("Cdag", 2)}),
      });

  // Spinful electrons.
  check_transpose_conj_hc_matrices(
      Electron(3),
      {
          OpSum(c, Op("Cdagup", 1)),
          OpSum(c, Op("Cdagdn", 2)),
          OpSum(c, Op("Cup", 0)),
          OpSum(c, Op("Cdn", 1)),
          OpSum(c, Op("Hop", {0, 2})),
          OpSum(c, Op("HopAsym", {0, 2})),
          OpSum(c, Op("Hopup", {0, 1})),
          OpSum(c, Op("Hopdn", {1, 2})),
          OpSum(c, Monomial{Op("Cdagup", 0), Op("Cdagdn", 2)}),
      });

  // tJ. HopAsym is not dispatched by the tJ kernels, so only Hop, Hopup and
  // Hopdn are checked here.
  check_transpose_conj_hc_matrices(
      tJ(3),
      {
          OpSum(c, Op("Cdagup", 1)),
          OpSum(c, Op("Cdagdn", 2)),
          OpSum(c, Op("Cup", 0)),
          OpSum(c, Op("Cdn", 1)),
          OpSum(c, Op("Hop", {0, 2})),
          OpSum(c, Op("Hopup", {0, 1})),
          OpSum(c, Op("Hopdn", {1, 2})),
          OpSum(c, Monomial{Op("Cdagup", 0), Op("Cdagdn", 2)}),
      });

  // --- Signs and partners, independent of the block ---
  // Cdag and C are transposes of each other, with sign +1
  REQUIRE(transpose(OpSum(Op("Cdag", 1))) == OpSum(Op("C", 1)));
  REQUIRE(transpose(OpSum(Op("Cdagup", 1))) == OpSum(Op("Cup", 1)));
  REQUIRE(transpose(OpSum(Op("Cdagdn", 2))) == OpSum(Op("Cdn", 2)));
  REQUIRE(transpose(OpSum(Op("Cup", 0))) == OpSum(Op("Cdagup", 0)));
  REQUIRE(transpose(OpSum(Op("Cdn", 1))) == OpSum(Op("Cdagdn", 1)));

  // Hop and Hopup are symmetric; HopAsym and HopupAsym are antisymmetric
  REQUIRE(transpose(OpSum(Op("Hop", {0, 2}))) == OpSum(Op("Hop", {0, 2})));
  REQUIRE(transpose(OpSum(Op("Hopup", {0, 1}))) == OpSum(Op("Hopup", {0, 1})));
  REQUIRE(transpose(OpSum(Op("HopAsym", {0, 2}))) ==
          OpSum(-1.0, Op("HopAsym", {0, 2})));
  REQUIRE(transpose(OpSum(Op("HopupAsym", {0, 1}))) ==
          OpSum(-1.0, Op("HopupAsym", {0, 1})));

  // Monomial: the transpose reverses the product of creation operators
  Monomial m{Op("Cdagup", 0), Op("Cdagdn", 2)};
  REQUIRE(transpose(m) == OpSum(Monomial{Op("Cdn", 2), Op("Cup", 0)}));

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}
