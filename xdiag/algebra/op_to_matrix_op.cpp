// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "op_to_matrix_op.hpp"

#include <vector>

#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::operators {

Op op_to_matrix_op(Op const &op) try {
  static const arma::mat sp = {{0., 1.}, {0., 0.}};
  static const arma::mat sm = {{0., 0.}, {1., 0.}};
  static const arma::mat sz = {{0.5, 0.}, {0., -0.5}};

  std::string const &type = op.type();

  if (type == "Matrix") {
    return op;
  }
  if (type == "S+") {
    return Op("Matrix", op.sites(), sp);
  }
  if (type == "S-") {
    return Op("Matrix", op.sites(), sm);
  }
  if (type == "Sz") {
    return Op("Matrix", op.sites(), sz);
  }

  // Sx{i} -> [[0,0.5],[0.5,0]]
  if (type == "Sx") {
    arma::mat m = {{0.0, 0.5}, {0.5, 0.0}};
    return Op("Matrix", op.sites(), m);
  }

  // Sy{i} -> [[0,-i/2],[i/2,0]]
  if (type == "Sy") {
    complex I(0.0, 1.0);
    arma::cx_mat m = {{0.0 + 0.0 * I, -I * 0.5}, {I * 0.5, 0.0 + 0.0 * I}};
    return Op("Matrix", op.sites(), m);
  }

  // Two-site compound ops. When i == j (same site) the matrix is a plain 2×2
  // product; otherwise it is a 4×4 Kronecker product. Convention for different
  // sites:  kron(outer=site1, inner=site0) — matches combine_matrix_ops.

  // SdotS{i,j} = Sz_i Sz_j + 0.5*S+_i S-_j + 0.5*S-_i S+_j
  if (type == "SdotS") {
    int64_t i = op.sites()[0], j = op.sites()[1];
    if (i == j) {
      arma::mat m = arma::mat(sz * sz) + 0.5 * arma::mat(sp * sm) +
                    0.5 * arma::mat(sm * sp);
      return Op("Matrix", std::vector<int64_t>{i}, m);
    }
    arma::mat m = arma::mat(arma::kron(sz, sz)) +
                  0.5 * arma::mat(arma::kron(sm, sp)) +
                  0.5 * arma::mat(arma::kron(sp, sm));
    return Op("Matrix", op.sites(), m);
  }

  if (type == "SzSz") {
    int64_t i = op.sites()[0], j = op.sites()[1];
    if (i == j) {
      return Op("Matrix", std::vector<int64_t>{i}, arma::mat(sz * sz));
    }
    return Op("Matrix", op.sites(), arma::mat(arma::kron(sz, sz)));
  }
  if (type == "Exchange") {
    int64_t i = op.sites()[0], j = op.sites()[1];
    if (i == j) {
      arma::mat m = 0.5 * arma::mat(sp * sm) + 0.5 * arma::mat(sm * sp);
      return Op("Matrix", std::vector<int64_t>{i}, m);
    }
    arma::mat m = 0.5 * arma::mat(arma::kron(sm, sp)) +
                  0.5 * arma::mat(arma::kron(sp, sm));
    return Op("Matrix", op.sites(), m);
  }

  // Three-site ScalarChirality → 8×8 complex matrix.
  // S_i·(S_j×S_k) = (i/2)*[ S+_i S-_j Sz_k - S-_i S+_j Sz_k
  //                        + Sz_i S+_j S-_k - Sz_i S-_j S+_k
  //                        + S-_i Sz_j S+_k - S+_i Sz_j S-_k ]
  // Each product A_i B_j C_k maps to kron(C_k, kron(B_j, A_i)).
  if (type == "ScalarChirality") {
    complex I(0., 1.);
    arma::cx_mat sp_cx(sp, arma::zeros<arma::mat>(2, 2));
    arma::cx_mat sm_cx(sm, arma::zeros<arma::mat>(2, 2));
    arma::cx_mat sz_cx(sz, arma::zeros<arma::mat>(2, 2));
    auto t = [](arma::cx_mat const &a, arma::cx_mat const &b,
                arma::cx_mat const &c) {
      return arma::cx_mat(arma::kron(a, arma::cx_mat(arma::kron(b, c))));
    };
    arma::cx_mat m = (I * 0.5) * t(sz_cx, sm_cx, sp_cx)     // S+_i S-_j Sz_k
                     + (-I * 0.5) * t(sz_cx, sp_cx, sm_cx)  // S-_i S+_j Sz_k
                     + (I * 0.5) * t(sm_cx, sp_cx, sz_cx)   // Sz_i S+_j S-_k
                     + (-I * 0.5) * t(sp_cx, sm_cx, sz_cx)  // Sz_i S-_j S+_k
                     + (I * 0.5) * t(sp_cx, sz_cx, sm_cx)   // S-_i Sz_j S+_k
                     + (-I * 0.5) * t(sm_cx, sz_cx, sp_cx); // S+_i Sz_j S-_k
    return Op("Matrix", op.sites(), m);
  }

  XDIAG_THROW(
      fmt::format("Cannot convert Op of type \"{}\" to a Matrix op", type));
} XDIAG_CATCH

} // namespace xdiag::operators
