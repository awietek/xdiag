//
// Created by Luke Staszewski on 13.06.23.
//
#pragma once
#include <utility>
#include "hydra/all.h"

template<typename mult, typename V>
inline std::pair<arma::cx_mat, arma::cx_mat> arnoldi_iteration(mult const& H,
                                                               V const& q0,
                                                               int n = 80,
                                                               double eps = 1e-12)
                                                               {
    /* following from implementation on wikipedia: https://en.wikipedia.org/wiki/Arnoldi_iteration
     * with edits based on Gram-Schmidt re-orthogonalisation: https://core.ac.uk/download/pdf/82406974.pdf
     * Parameters:
     * -----------
     * H: hamiltonian, (m, m) given as lambda function (H(v, w) applies H to v, stored as w)
     * q0: initial vector for Arnoldi process (m,)
     * n: number of iterations - n.b. should stop if n> m-1
     * eps: deflation tolerance for stopping the process i.e. if beta < eps
     *
     * Returns:
     * --------
     * Q: matrix of Arnoldi vectors (m, n+1)
     * h: hamiltonian in the basis of the Arnoldi vectors (n+1, n) (n.b. may want to discard final row)
     * TODO: add convergence check
     * TODO: write to disk
     */

    using namespace arma;
    using namespace hydra;

    int m = q0.n_elem;

    cx_mat h(n+1, n, fill::zeros);
    cx_mat Q(m, n+1);

    // normalise q0 and add to Q
    Q.col(0) = q0/norm(q0);

    // Arnoldi iterations
    for (int k = 1; k < n+1; ++k){
        cx_vec v;
        H(Q.col(k-1), v);
        auto v_copy = v;
        for (int j = 0; j < k; ++j) {
            h(j, k - 1) = cdot(Q.col(j), v);
            v -= h(j, k - 1) * Q.col(j);
        }
        for (int j = 0; j < k; ++j) {
            auto temp = cdot(Q.col(j), v);
            v -= temp * Q.col(j);
        }

        h(k, k-1) = norm(v);
        // check for breakdown
        if (abs(h(k, k-1)) < eps){
            return {Q(span(0, m-1), span(0, k-1)),
                    h(span(0, k), span(0, k-1))};
        }
        Q.col(k) = v/ h(k, k-1);
    } // end of Arnoldi iterations

    return {Q(span(0, m-1), span(0, n-1)),
            h};
}


template<typename mult, typename V>
inline std::pair<arma::cx_vec, arma::cx_mat> arnoldi(mult const& H,
                                                  V const& q0,
                                                  int n, double eps = 1e-12,
                                                  bool is_hermitian = false)
{
    /* Arnoldi scheme for eigenvector decomposition
     * Parameters
     * ----------
     * H: matrix to be diagonalised
     * q0: starting vector for Arnoldi scheme (see krylov space)
     * n: iterations for Arnoldi
     * eps: breakdown tolerance of Arnoldi
     * is_hermitian: if H is hermitian option for eig_sym in algorithm for improved performance
     * Returns
     * -------
     * eigvals: vec
     * eigvecs: cx_mat
     */

    using namespace hydra;
    using namespace arma;
    auto arnoldi_out = arnoldi_iteration(H, q0, n, eps);
    auto Q = arnoldi_out.first;
    auto h = arnoldi_out.second;

    // krylov space dimension
    Log.out(1, "h_rows: {}, h_cols: {}", h.n_rows, h.n_cols);
    Log.out(1, "Q_rows: {}, Q_cols: {}", Q.n_rows, Q.n_cols);
    int M = Q.n_cols; // krylov space dim
    uword m = q0.n_elem;
    h = h(span(0, M-1), span(0, M-1));
    cx_vec ritz_eigvals;
    cx_mat ritz_eigvecs;
    if (is_hermitian){
        vec eigs_real;
        eig_sym(eigs_real, ritz_eigvecs, h);
        ritz_eigvals = conv_to<cx_vec>::from(eigs_real);
        ritz_eigvecs = Q*ritz_eigvecs;
        return {ritz_eigvals, ritz_eigvecs};
    }
    eig_gen(ritz_eigvals, ritz_eigvecs, h);

    // now finding the eigenvectors in the original basis
    ritz_eigvecs = Q*ritz_eigvecs;
    return {ritz_eigvals, ritz_eigvecs};
}

