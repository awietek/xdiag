//
// Created by Luke Staszewski on 13.06.23.
//
#pragma once
#include <hydra/extern/armadillo/armadillo>
#include <hydra/extern/fmt/format.h>
#include <hydra/utils/logger.h>
#include <hydra/utils/timing.h>
#include <string>
#include <utility>

namespace hydra {

template <typename mult, typename coeff_t>
inline bool arnoldi_step(mult const &H, arma::Mat<coeff_t> &h,
                         std::string path_to_Avecs, int k, double eps = 1e-12) {
  /*
   * Returns
   * -------
   * bool: true if breakdown has occured
   */
  using namespace arma;
  using vec_t = Col<coeff_t>;

  vec_t v_prev;
  v_prev.load(fmt::format("{}/Arnoldi_{}.arm", path_to_Avecs, k - 1));
  vec_t v(v_prev);
  auto t0 = rightnow();
  H(v_prev, v);
  timing(rightnow(), t0, "Arnoldi MVM", 2);

  t0 = rightnow();
  for (int j = 0; j < k; ++j) {
    v_prev.load(fmt::format("{}/Arnoldi_{}.arm", path_to_Avecs, j));
    h(j, k - 1) = cdot(v_prev, v);
    v -= h(j, k - 1) * v_prev;
  }
  for (int j = 0; j < k; ++j) {
    v_prev.load(fmt::format("{}/Arnoldi_{}.arm", path_to_Avecs, j));
    auto temp = cdot(v_prev, v);
    v -= temp * v_prev;
  }
  timing(rightnow(), t0, "Arnoldi ortho", 2);
  
  h(k, k - 1) = norm(v);
  // check for breakdown
  if (abs(h(k, k - 1)) < eps) {
    return true;
  }
  v /= h(k, k - 1);
  v.save(fmt::format("{}/Arnoldi_{}.arm", path_to_Avecs, k));

  return false;
}

template <typename mult, typename coeff_t>
inline arma::Mat<coeff_t>
arnoldi_iteration(mult const &H, arma::Col<coeff_t> const &q0,
                  std::string path_to_Avecs, int n = 80, double eps = 1e-12) {
  /* following from implementation on wikipedia:
   * https://en.wikipedia.org/wiki/Arnoldi_iteration with edits based on
   * Gram-Schmidt re-orthogonalisation:
   * https://core.ac.uk/download/pdf/82406974.pdf Parameters:
   * -----------
   * H: hamiltonian, (m, m) given as lambda function (H(v, w) applies H to v,
   * stored as w) q0: initial vector for Arnoldi process (m,) path_to_Avecs:
   * where the Arnoldi vectors will be stored (abs.) n: number of iterations -
   * n.b. should stop if n> m-1 eps: deflation tolerance for stopping the
   * process i.e. if beta < eps
   *
   * Returns:
   * --------
   * h: hamiltonian in the basis of the Arnoldi vectors (n+1, n) (n.b. may want
   * to discard final row)
   * TODO: add convergence check
   */
  using namespace arma;
  using mat_t = Mat<coeff_t>;
  using vec_t = Col<coeff_t>;

  mat_t h(n + 1, n, fill::zeros);

  // normalise q0 and add to Q
  double norm_0 = norm(q0);
  if (norm_0 == 0) {
    Log.err("initial vector norm is zero");
    exit(1);
  }
  vec_t v = q0 / norm(q0);
  v.save(fmt::format("{}/Arnoldi_{}.arm", path_to_Avecs, 0));

  bool breakdown = false;
  // Arnoldi iterations
  for (int k = 1; k < n + 1; ++k) {
    breakdown = arnoldi_step(H, h, path_to_Avecs, k, eps);
    if (breakdown) {
      return h(span(0, k), span(0, k - 1));
    }
  } // end of Arnoldi iterations

  return h;
}

template <typename mult, typename coeff_t>
inline std::pair<arma::Col<coeff_t>, arma::Mat<coeff_t>>
arnoldi(mult const &H, arma::Col<coeff_t> const &q0, std::string path_to_vecs,
        int n, double eps = 1e-12, bool is_hermitian = false,
        bool build_ritz_vecs = true) {
  /* Arnoldi scheme for eigenvector decomposition
   * Parameters
   * ----------
   * H: matrix to be diagonalised
   * q0: starting vector for Arnoldi scheme (see krylov space)
   * path_to_vecs: for storage of eigenvectors and Arnoldi vectors
   * n: iterations for Arnoldi
   * eps: breakdown tolerance of Arnoldi
   * is_hermitian: if H is hermitian option for eig_sym in algorithm for
   * improved performance Returns
   * -------
   * eigvals: vec
   * n.b. eigvecs will be stored along with Arnoldi vecs in path_to_vecs
   * TODO: template out eig vals being real or complex
   */

  using namespace arma;
  using mat_t = Mat<coeff_t>;
  using vec_t = Col<coeff_t>;

  auto h = arnoldi_iteration(H, q0, path_to_vecs, n, eps);

  // krylov space dimension
  int M = h.n_cols; // krylov space dim
  uword m = q0.n_elem;
  h = h(span(0, M - 1), span(0, M - 1));
  vec_t ritz_eigvals;

  if (build_ritz_vecs) {
    mat_t ritz_eigvecs_Abasis; // Arnoldi basis
    if (is_hermitian) {
      vec eigs_real;
      eig_sym(eigs_real, ritz_eigvecs_Abasis, h);
      ritz_eigvals = conv_to<vec_t>::from(eigs_real);
      // need to build up each eigenvector
      for (int i = 0; i < M; ++i) {
        vec_t eig_vec_i(m);
        vec_t arnoldi_k(m);
        for (int k = 0; k < M; ++k) {
          arnoldi_k.load(fmt::format("{}/Arnoldi_{}.arm", path_to_vecs, k));
          eig_vec_i += arnoldi_k * ritz_eigvecs_Abasis(k, i);
        }
        eig_vec_i.save(fmt::format("{}/Ritz_{}.arm", path_to_vecs, i));
      }
      return {ritz_eigvals, h};
    }

    eig_gen(ritz_eigvals, ritz_eigvecs_Abasis, h);
    // need to build up each eigenvector
    for (int i = 0; i < M; ++i) {
      vec_t eig_vec_i(m);
      vec_t arnoldi_k(m);
      for (int k = 0; k < M; ++k) {
        arnoldi_k.load(fmt::format("{}/Arnoldi_{}.arm", path_to_vecs, k));
        eig_vec_i += arnoldi_k * ritz_eigvecs_Abasis(k, i);
      }
      eig_vec_i.save(fmt::format("{}/Ritz_{}.arm", path_to_vecs, i));
    }
  }

  return {ritz_eigvals, h};
}

template <typename coeff_t>
arma::Mat<coeff_t> read_arnoldi_vectors(std::string path_to_vecs, int n = 0);
arma::mat read_arnoldi_vectors_real(std::string path_to_vecs, int n = 0);
arma::cx_mat read_arnoldi_vectors_cplx(std::string path_to_vecs, int n = 0);

template <typename coeff_t>
arma::Mat<coeff_t> read_ritz_vectors(std::string path_to_vecs, int n = 0);
arma::mat read_ritz_vectors_real(std::string path_to_vecs, int n = 0);
arma::cx_mat read_ritz_vectors_cplx(std::string path_to_vecs, int n = 0);

} // namespace hydra
