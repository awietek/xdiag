// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "op_to_matrix_op.hpp"

#include <cmath>
#include <vector>

#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::algebra {

namespace {

// Spin-S elementary matrices for local dimension d = 2S+1. Basis index
// i = 0..d-1 maps to magnetic quantum number m_i = -S + i, so index 0 is the
// lowest-weight state m = -S. This matches the occupation/bit convention used
// by the basis and the single-site kernels (a set bit is "up", i.e. higher m),
// so that the Matrix path for products agrees with term_sz / term_spsm. For
// d = 2 it gives Sz = diag(-1/2, 1/2), S+ = [[0,0],[1,0]], S- = [[0,1],[0,0]].
struct SpinMatrices {
  arma::mat sp;
  arma::mat sm;
  arma::mat sz;
};

SpinMatrices spin_matrices(int64_t d) {
  double S = 0.5 * static_cast<double>(d - 1);
  arma::mat sp(d, d, arma::fill::zeros);
  arma::mat sm(d, d, arma::fill::zeros);
  arma::mat sz(d, d, arma::fill::zeros);

  for (int64_t i = 0; i < d; ++i) {
    sz(i, i) = -S + static_cast<double>(i);
  }
  // S+ |m> = sqrt(S(S+1) - m(m+1)) |m+1>; |m+1> sits at index i+1.
  for (int64_t i = 0; i < d - 1; ++i) {
    double m = -S + static_cast<double>(i);
    sp(i + 1, i) = std::sqrt(S * (S + 1.0) - m * (m + 1.0));
  }
  // S- |m> = sqrt(S(S+1) - m(m-1)) |m-1>; |m-1> sits at index i-1.
  for (int64_t i = 1; i < d; ++i) {
    double m = -S + static_cast<double>(i);
    sm(i - 1, i) = std::sqrt(S * (S + 1.0) - m * (m - 1.0));
  }
  return {sp, sm, sz};
}

} // namespace

Op op_to_matrix_op(Op const &op, int64_t d) try {
  std::string const &type = op.type();

  if (type == "Matrix") {
    return op;
  }

  // --- Bosonic operators on a Fock space truncated to occupations 0..d-1 -----
  // Basis index i = occupation i.
  //   A    |n> = sqrt(n)   |n-1>   (annihilation)
  //   Adag |n> = sqrt(n+1) |n+1>   (creation, truncated at n = d-1)
  //   N    |n> = n |n>             (number)
  if (type == "A" || type == "Adag" || type == "N") {
    arma::mat m(d, d, arma::fill::zeros);
    if (type == "N") {
      for (int64_t n = 0; n < d; ++n) {
        m(n, n) = static_cast<double>(n);
      }
    } else if (type == "A") {
      for (int64_t n = 1; n < d; ++n) {
        m(n - 1, n) = std::sqrt(static_cast<double>(n));
      }
    } else { // Adag
      for (int64_t n = 0; n < d - 1; ++n) {
        m(n + 1, n) = std::sqrt(static_cast<double>(n + 1));
      }
    }
    return Op("Matrix", op.sites(), m);
  }

  // --- Spin-S operators (S = (d-1)/2) ----------------------------------------
  SpinMatrices spin = spin_matrices(d);
  arma::mat const &sp = spin.sp;
  arma::mat const &sm = spin.sm;
  arma::mat const &sz = spin.sz;

  if (type == "S+") {
    return Op("Matrix", op.sites(), sp);
  }
  if (type == "S-") {
    return Op("Matrix", op.sites(), sm);
  }
  if (type == "Sz") {
    return Op("Matrix", op.sites(), sz);
  }

  // Sx = (S+ + S-) / 2
  if (type == "Sx") {
    arma::mat m = 0.5 * (sp + sm);
    return Op("Matrix", op.sites(), m);
  }

  // Sy = (S+ - S-) / (2i) = -(i/2) (S+ - S-)
  if (type == "Sy") {
    complex I(0.0, 1.0);
    arma::cx_mat m =
        arma::cx_mat(arma::zeros<arma::mat>(d, d), -0.5 * (sp - sm));
    return Op("Matrix", op.sites(), m);
  }

  // Two-site compound ops. When i == j (same site) the matrix is a plain d×d
  // product; otherwise it is a d²×d² Kronecker product. Convention for
  // different sites:  kron(outer=site1, inner=site0) — matches
  // combine_matrix_ops.

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
  if (type == "ExchangeAsym") {
    int64_t i = op.sites()[0], j = op.sites()[1];
    if (i == j) {
      arma::mat m = 0.5 * arma::mat(sp * sm) - 0.5 * arma::mat(sm * sp);
      return Op("Matrix", std::vector<int64_t>{i}, m);
    }
    arma::mat m = 0.5 * arma::mat(arma::kron(sm, sp)) -
                  0.5 * arma::mat(arma::kron(sp, sm));
    return Op("Matrix", op.sites(), m);
  }

  // Three-site ScalarChirality → d³×d³ complex matrix.
  // S_i·(S_j×S_k) = (i/2)*[ S+_i S-_j Sz_k - S-_i S+_j Sz_k
  //                        + Sz_i S+_j S-_k - Sz_i S-_j S+_k
  //                        + S-_i Sz_j S+_k - S+_i Sz_j S-_k ]
  // Each product A_i B_j C_k maps to kron(C_k, kron(B_j, A_i)).
  if (type == "ScalarChirality") {
    complex I(0., 1.);
    arma::mat z = arma::zeros<arma::mat>(d, d);
    arma::cx_mat sp_cx(sp, z);
    arma::cx_mat sm_cx(sm, z);
    arma::cx_mat sz_cx(sz, z);
    auto t = [](arma::cx_mat const &a, arma::cx_mat const &b,
                arma::cx_mat const &c) {
      return arma::cx_mat(arma::kron(a, arma::cx_mat(arma::kron(b, c))));
    };
    // Realizes ScalarChirality = S_i . (S_j x S_k); matches the
    // term_scalar_chirality kernel so the kernel and Matrix path agree.
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
}
XDIAG_CATCH

} // namespace xdiag::algebra
