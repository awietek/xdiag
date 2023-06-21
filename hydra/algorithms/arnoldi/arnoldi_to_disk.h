//
// Created by Luke Staszewski on 13.06.23.
//
#pragma once
#include <utility>
#include "hydra/all.h"

template<typename mult>
inline bool arnoldi_step(mult const& H, arma::cx_mat & h,
                         std::string path_to_Avecs,
                         int k, double eps = 1e-12)
                         {
    /*
     * Returns
     * -------
     * bool: true if breakdown has occured
     */
    using namespace arma;
    using namespace hydra;

    cx_vec v_prev;
    v_prev.load(fmt::format("{}/Arnoldi_{}.bin", path_to_Avecs, k-1));
    cx_vec v(v_prev);
    H(v_prev, v);
    for (int j = 0; j < k; ++j) {
        v_prev.load(fmt::format("{}/Arnoldi_{}.bin", path_to_Avecs, j));
        h(j, k - 1) = cdot(v_prev, v);
        v -= h(j, k - 1) * v_prev;
    }
    for (int j = 0; j < k; ++j) {
        v_prev.load(fmt::format("{}/Arnoldi_{}.bin", path_to_Avecs, j));
        auto temp = cdot(v_prev, v);
        v -= temp * v_prev;
    }

    h(k, k-1) = norm(v);
    // check for breakdown
    if (abs(h(k, k-1)) < eps){
        return true;
    }
    v /= h(k, k-1);
    v.save(fmt::format("{}/Arnoldi_{}.bin", path_to_Avecs, k));

    return false;
}


template<typename mult>
inline arma::cx_mat arnoldi_iteration(mult const& H, arma::cx_vec const& q0, std::string path_to_Avecs, int n = 80, double eps = 1e-12)
{
    /* following from implementation on wikipedia: https://en.wikipedia.org/wiki/Arnoldi_iteration
     * with edits based on Gram-Schmidt re-orthogonalisation: https://core.ac.uk/download/pdf/82406974.pdf
     * Parameters:
     * -----------
     * H: hamiltonian, (m, m) given as lambda function (H(v, w) applies H to v, stored as w)
     * q0: initial vector for Arnoldi process (m,)
     * path_to_Avecs: where the Arnoldi vectors will be stored (abs.)
     * n: number of iterations - n.b. should stop if n> m-1
     * eps: deflation tolerance for stopping the process i.e. if beta < eps
     *
     * Returns:
     * --------
     * h: hamiltonian in the basis of the Arnoldi vectors (n+1, n) (n.b. may want to discard final row)
     * TODO: add convergence check
     */

    using namespace arma;
    using namespace hydra;

    cx_mat h(n+1, n, fill::zeros);

    // normalise q0 and add to Q
    double norm_0 = norm(q0);
    if (norm_0 == 0){
        Log.err("initial vector norm is zero");
        exit(1);
    }
    cx_vec v = q0/norm(q0);
    v.save(fmt::format("{}/Arnoldi_{}.bin", path_to_Avecs, 0));

    bool breakdown = false;
    // Arnoldi iterations
    for (int k = 1; k < n+1; ++k){
        breakdown = arnoldi_step(H, h, path_to_Avecs, k, eps);
        if (breakdown){
            return h(span(0, k), span(0, k-1));
        }
    } // end of Arnoldi iterations

    return h;
}


template<typename mult>
inline arma::cx_vec arnoldi(mult const& H,
                            arma::cx_vec const& q0,
                            std::string path_to_vecs,
                            int n, double eps = 1e-12,
                            bool is_hermitian = false)
{
    /* Arnoldi scheme for eigenvector decomposition
     * Parameters
     * ----------
     * H: matrix to be diagonalised
     * q0: starting vector for Arnoldi scheme (see krylov space)
     * path_to_vecs: for storage of eigenvectors and Arnoldi vectors
     * n: iterations for Arnoldi
     * eps: breakdown tolerance of Arnoldi
     * is_hermitian: if H is hermitian option for eig_sym in algorithm for improved performance
     * Returns
     * -------
     * eigvals: vec
     * n.b. eigvecs will be stored along with Arnoldi vecs in path_to_vecs
     * TODO: template out eig vals being real or complex
     */

    using namespace hydra;
    using namespace arma;
    auto h = arnoldi_iteration(H, q0, path_to_vecs, n, eps);

    // krylov space dimension
    int M = h.n_cols; // krylov space dim
    uword m = q0.n_elem;
    h = h(span(0, M-1), span(0, M-1));
    cx_vec ritz_eigvals;
    cx_mat ritz_eigvecs_Abasis; // Arnoldi basis
    if (is_hermitian){
        vec eigs_real;
        eig_sym(eigs_real, ritz_eigvecs_Abasis, h);
        ritz_eigvals = conv_to<cx_vec>::from(eigs_real);
        // need to build up each eigenvector
        for (int i = 0; i < M; ++i) {
            cx_vec eig_vec_i(m);
            cx_vec arnoldi_k(m);
            for (int k = 0; k < M; ++k) {
                arnoldi_k.load(fmt::format("{}/Arnoldi_{}.bin", path_to_vecs, k));
                eig_vec_i += arnoldi_k * ritz_eigvecs_Abasis(k, i);
            }
            eig_vec_i.save(fmt::format("{}/Ritz_{}.bin", path_to_vecs, i));
        }
        return ritz_eigvals;
    }

    eig_gen(ritz_eigvals, ritz_eigvecs_Abasis, h);
    // need to build up each eigenvector
    for (int i = 0; i < M; ++i) {
        cx_vec eig_vec_i(m);
        cx_vec arnoldi_k(m);
        for (int k = 0; k < M; ++k) {
            arnoldi_k.load(fmt::format("{}/Arnoldi_{}.bin", path_to_vecs, k));
            eig_vec_i += arnoldi_k * ritz_eigvecs_Abasis(k, i);
        }
        eig_vec_i.save(fmt::format("{}/Ritz_{}.bin", path_to_vecs, i));
    }

    return ritz_eigvals;
}

