//
// Created by Luke Staszewski on 19.06.23.
//
#include "../../catch.hpp"
#include <iostream>

#include "../../blocks/electron/testcases_electron.hpp"

#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/arnoldi/arnoldi.hpp>
#include <xdiag/algorithms/arnoldi/arnoldi_to_disk.hpp>
#include <xdiag/algorithms/time_evolution/pade_matrix_exponential.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/random_state.hpp>

// TODO: write tests for arnoldi_to_disk
bool check_basis_orthonormality(arma::cx_mat const &Q, double tol = 1e-12) {
  using namespace arma;
  int M = Q.n_cols;
  cx_mat Q_eye = trans(Q) * Q;
  cx_mat diff = Q_eye - eye(M, M);
  if (norm(diff) < tol) {
    return true;
  }
  return false;
}

template <typename mult, typename B>
bool check_eigenvector_equation(mult const &H, arma::cx_mat eigvecs, B eigvals,
                                double tol = 1e-12) {
  using namespace arma;
  uword M = eigvecs.n_cols;
  for (int i = 0; i < M; ++i) {
    arma::cx_vec w;
    H(eigvecs.col(i), w);
    double diff = arma::norm(w - eigvals(i) * eigvecs.col(i));
    if (diff > tol) {
      xdiag::Log.warn("eigenvector {} doesn't satisfy eigvec eq., diff = {}", i,
                      diff);
      return false;
    }
  }
  return true;
}

TEST_CASE("ritz_vecs_arnoldi", "[arnoldi]") {
  using namespace xdiag;
  using namespace arma;

  Log("testing ritz_vecs_arnoldi");

  // test getting eigenvalues/ eigenvectors for (U + U_dagg) and H
  int n_sites = 6;
  auto block = Spinhalf(n_sites);
  OpSum ops;

  for (int i = 0; i < n_sites; ++i) {
    ops += "J1" * Op("SdotS", {i, (i + 1) % n_sites});
    ops += "J2" * Op("SdotS", {i, (i + 2) % n_sites});
  }
  ops["J1"] = 1;
  ops["J2"] = .5;

  auto state_0 = rand(block, false);
  cx_mat A = complex(0, -1) * matrixC(ops, block);
  cx_mat U = expm(A);
  cx_mat T = trans(U) + U;

  auto apply_T = [&T](cx_vec const &v, cx_vec &w) { w = T * v; };

  {
    // getting the eigenvalues through arnoldi
    auto arnoldi_out = arnoldi(apply_T, state_0.vectorC(), 70, 1e-12, false);
    auto ritz_vals = arnoldi_out.first;
    auto ritz_vecs = arnoldi_out.second;

    REQUIRE(check_eigenvector_equation(apply_T, ritz_vecs, ritz_vals, 1e-12));

    auto arnoldi_out_sym = arnoldi(apply_T, state_0.vectorC(), 70, 1e-12, true);
    auto ritz_vals_sym = arnoldi_out_sym.first;
    auto ritz_vecs_sym = arnoldi_out_sym.second;

    REQUIRE(check_eigenvector_equation(apply_T, ritz_vecs_sym, ritz_vals_sym,
                                       1e-12));
    REQUIRE(check_basis_orthonormality(ritz_vecs_sym));
    REQUIRE(check_eigenvector_equation(apply_T, ritz_vecs_sym, ritz_vals_sym,
                                       1e-12));
  }

  {
    std::string dumpdir = XDIAG_DIRECTORY "/misc/dump";
    // getting the eigenvalues through arnoldi
    auto [ritz_vals, h] =
        arnoldi(apply_T, state_0.vectorC(), dumpdir, 70, 1e-12, false);

    auto arno_vecs = read_arnoldi_vectors_cplx(dumpdir);

    auto ritz_vecs = read_ritz_vectors_cplx(dumpdir);
    REQUIRE(check_eigenvector_equation(apply_T, ritz_vecs, ritz_vals, 1e-12));
    auto [ritz_vals_sym, h_sym] =
        arnoldi(apply_T, state_0.vectorC(), dumpdir, 70, 1e-12, true);
    auto ritz_vecs_sym = read_ritz_vectors_cplx(dumpdir);
    auto arno_vecs_sym = read_arnoldi_vectors_cplx(dumpdir);

    REQUIRE(check_eigenvector_equation(apply_T, ritz_vecs_sym, ritz_vals_sym,
                                       1e-12));
    REQUIRE(check_basis_orthonormality(arno_vecs_sym));
    REQUIRE(check_basis_orthonormality(ritz_vecs_sym));
    REQUIRE(check_eigenvector_equation(apply_T, ritz_vecs_sym, ritz_vals_sym,
                                       1e-12));
  }
}
