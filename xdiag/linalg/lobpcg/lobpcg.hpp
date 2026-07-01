// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/math/ipow.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

// Locally Optimal Block Preconditioned Conjugate Gradient (LOBPCG), specialized
// to the *standard* Hermitian eigenproblem A X = X Lambda (no mass matrix B, no
// preconditioner, no constraints). Computes the "neigs" algebraically smallest
// eigenvalues. A block of size "neigs + guard" is iterated so that a degenerate
// multiplet straddling the neigs-th eigenvalue is captured with the correct
// multiplicity (the returned lowest-neigs Ritz values are exact eigenvalues at
// convergence, unlike single-vector Lanczos).
//
// Reference: A. V. Knyazev, SIAM J. Sci. Comput. 23 (2001) 517.
// Structure follows scipy.sparse.linalg.lobpcg with B = M = Y = None.

namespace xdiag::linalg {

struct lobpcg_result_t {
  arma::vec eigenvalues;    // block-size Ritz values, ascending
  arma::vec residual_norms; // 2-norm of A x - lambda x, per column
  int64_t niterations;
  std::string criterion; // "converged" | "maxiterations" | "failed"
  // Per-iteration history (one row per iteration, blocksize columns). Cheap
  // (small dense) and useful for convergence/degeneracy diagnostics.
  arma::mat eigenvalue_history;
  arma::mat residual_norms_history;
};

// Small dense generalized Hermitian eigenproblem  A z = lambda B z, with B SPD.
// Reduces to a standard problem via the Cholesky factor of B (see
// https://www.netlib.org/lapack/lug/node54.html):
//   B = L L^H,  C = L^-1 A L^-H,  C w = lambda w,  z = L^-H w.
// The resulting z are B-orthonormal (z^H B z = I), which is what the
// Rayleigh-Ritz step requires. Returns false (never throws) on numerical
// failure so the caller can restart with a smaller trial subspace.
template <typename coeff_t>
bool eig_sym_gen(arma::vec &eigvals, arma::Mat<coeff_t> &eigvecs,
                 arma::Mat<coeff_t> const &A, arma::Mat<coeff_t> const &B) {
  arma::Mat<coeff_t> L;
  if (!arma::chol(L, B, "lower")) { // B = L L^H, L lower triangular
    return false;
  }
  arma::Mat<coeff_t> Linv;
  if (!arma::inv(Linv, arma::trimatl(L))) {
    return false;
  }
  // C = L^-1 A L^-H (.t() is the conjugate transpose, so Linv.t() == L^-H)
  arma::Mat<coeff_t> C = Linv * A * Linv.t();
  C = (C + C.t()) / 2.0; // symmetrize away roundoff before the solver
  if (!arma::eig_sym(eigvals, eigvecs, C)) {
    return false;
  }
  // recover the generalized eigenvectors z = L^-H w by solving L^H z = w
  arma::Mat<coeff_t> Z;
  if (!arma::solve(Z, arma::trimatu(L.t()), eigvecs)) {
    return false;
  }
  eigvecs = std::move(Z);
  return true;
}

// In-place orthonormalization of a block via Cholesky-QR: given V, factor
// V^H V = R^H R and set V <- V R^-1 so that V^H V = I afterwards. Returns the
// factor R^-1 in "invR" so the caller can update an associated block, e.g.
// A V, by the same right-multiplication (A V) R^-1. Returns false if the
// Cholesky fails (V numerically rank-deficient) so the caller can restart.
template <typename coeff_t, class dot_f>
bool orthonormalize(arma::Mat<coeff_t> &V, arma::Mat<coeff_t> &invR,
                    dot_f dot) {
  arma::Mat<coeff_t> gram = dot(V, V);
  gram = (gram + gram.t()) / 2.0;
  arma::Mat<coeff_t> R;
  if (!arma::chol(R, gram)) { // gram = R^H R, R upper triangular
    return false;
  }
  if (!arma::inv(invR, arma::trimatu(R))) {
    return false;
  }
  V = V * invR;
  return true;
}

// Assemble a symmetric 2x2 block matrix [[A, B], [B^H, D]].
template <typename coeff_t>
arma::Mat<coeff_t> sym_block_2x2(arma::Mat<coeff_t> const &A,
                                 arma::Mat<coeff_t> const &B,
                                 arma::Mat<coeff_t> const &D) {
  int64_t n1 = A.n_rows;
  int64_t n2 = D.n_rows;
  arma::Mat<coeff_t> M(n1 + n2, n1 + n2, arma::fill::zeros);
  using arma::span;
  M.submat(span(0, n1 - 1), span(0, n1 - 1)) = A;
  M.submat(span(0, n1 - 1), span(n1, n1 + n2 - 1)) = B;
  M.submat(span(n1, n1 + n2 - 1), span(0, n1 - 1)) = B.t();
  M.submat(span(n1, n1 + n2 - 1), span(n1, n1 + n2 - 1)) = D;
  return M;
}

// Assemble a symmetric 3x3 block matrix from the upper-triangular blocks
// [[A, B, C], [B^H, E, F], [C^H, F^H, I]].
template <typename coeff_t>
arma::Mat<coeff_t>
sym_block_3x3(arma::Mat<coeff_t> const &A, arma::Mat<coeff_t> const &B,
              arma::Mat<coeff_t> const &C, arma::Mat<coeff_t> const &E,
              arma::Mat<coeff_t> const &F, arma::Mat<coeff_t> const &I) {
  int64_t n1 = A.n_rows;
  int64_t n2 = E.n_rows;
  int64_t n3 = I.n_rows;
  int64_t n = n1 + n2 + n3;
  arma::Mat<coeff_t> M(n, n, arma::fill::zeros);
  using arma::span;
  M.submat(span(0, n1 - 1), span(0, n1 - 1)) = A;
  M.submat(span(0, n1 - 1), span(n1, n1 + n2 - 1)) = B;
  M.submat(span(0, n1 - 1), span(n1 + n2, n - 1)) = C;
  M.submat(span(n1, n1 + n2 - 1), span(0, n1 - 1)) = B.t();
  M.submat(span(n1, n1 + n2 - 1), span(n1, n1 + n2 - 1)) = E;
  M.submat(span(n1, n1 + n2 - 1), span(n1 + n2, n - 1)) = F;
  M.submat(span(n1 + n2, n - 1), span(0, n1 - 1)) = C.t();
  M.submat(span(n1 + n2, n - 1), span(n1, n1 + n2 - 1)) = F.t();
  M.submat(span(n1 + n2, n - 1), span(n1 + n2, n - 1)) = I;
  return M;
}

// 2-norm of each column of R, using the (possibly distributed) block dot
// product dot(R,R) whose diagonal holds the squared column norms.
template <typename coeff_t, class dot_f>
arma::vec column_norms(arma::Mat<coeff_t> const &R, dot_f dot) {
  arma::Mat<coeff_t> g = dot(R, R);
  // diagonal entries are the (real, non-negative) squared column norms; abs()
  // yields a real vector for both real and complex coeff_t.
  arma::Col<coeff_t> d = g.diag();
  return arma::sqrt(arma::abs(d));
}

// LOBPCG for the standard eigenproblem A X = X Lambda.
//
// multiplyA(V, W): compute W = A * V for n x m blocks V, W (the matrix-block
//                  operation; the only place the operator A enters).
// dot(V, W):       return V^H W as a small dense matrix. For serial blocks this
//                  is V.t()*W; for distributed blocks it must perform the
//                  cross-rank reduction (xdiag::math::matrix_dot).
// X:               n x blocksize initial guess (blocksize = neigs + guard),
//                  overwritten with the Ritz vectors on output.
// neigs:           number of lowest eigenvalues that must reach "tol".
//
// Large n x blocksize arrays held: X (in/out), AX, P, AP, R, AR, bestX, and two
// update temporaries -> one matrix-block MVM per iteration.
template <typename coeff_t, class multiply_f, class dot_f>
lobpcg_result_t lobpcg(multiply_f multiplyA, dot_f dot, arma::Mat<coeff_t> &X,
                       int64_t neigs, double tol = 1e-10,
                       int64_t max_iterations = 1000,
                       int64_t restart_control = 20) try {
  using mat_t = arma::Mat<coeff_t>;

  int64_t n = X.n_rows;
  int64_t blocksize = X.n_cols;
  if ((n == 0) || (blocksize == 0)) {
    XDIAG_THROW("Input initial block X has zero dimension");
  }
  if (neigs < 1) {
    XDIAG_THROW("Number of eigenvalues \"neigs\" must be >= 1");
  }
  if (neigs > blocksize) {
    XDIAG_THROW("Block size (columns of X) must be >= neigs");
  }
  // Note: n here is X.n_rows, which for a distributed block is the *local*
  // slice size, not the global dimension. The "problem too small for the block
  // size" reliability heuristic (dim >= 5*blocksize) therefore belongs to the
  // caller, which knows the global dimension -- it is not enforced here.
  if (tol <= 0.0) {
    tol = std::sqrt(std::numeric_limits<double>::epsilon()) * blocksize;
  }

  Log(1, "LOBPCG: matrix size {}, block size {}, neigs {}", n, blocksize,
      neigs);

  mat_t invR; // Cholesky factor returned by orthonormalize (reused)

  // Orthonormalize the initial block and compute the initial Ritz vectors.
  if (!orthonormalize(X, invR, dot)) {
    XDIAG_THROW("Linearly dependent initial approximations in LOBPCG");
  }
  mat_t AX(n, blocksize);
  multiplyA(X, AX);
  mat_t gramXAX = dot(X, AX);
  gramXAX = (gramXAX + gramXAX.t()) / 2.0;
  arma::vec lambda;
  mat_t coeff;
  if (!arma::eig_sym(lambda, coeff, gramXAX)) {
    XDIAG_THROW("Initial Rayleigh-Ritz eigensolve failed");
  }
  X = X * coeff;   // rotate to Ritz vectors (columns ascending in lambda)
  AX = AX * coeff; // keep A X consistent without a new MVM

  // Persistent auxiliary blocks (allocated once, reused each iteration).
  mat_t R(n, blocksize, arma::fill::zeros);  // residuals
  mat_t AR(n, blocksize, arma::fill::zeros); // A * (active residuals)
  mat_t P(n, blocksize, arma::fill::zeros);  // implicit conjugate directions
  mat_t AP(n, blocksize, arma::fill::zeros); // A * P
  mat_t bestX = X;                           // smallest-residual iterate
  arma::vec bestLambda = lambda;

  std::vector<arma::rowvec> lambdaHistory;
  std::vector<arma::rowvec> residualHistory;

  // Soft locking: once a column's residual drops below tol it is "locked" and
  // excluded from the search directions R, P (monotone). X keeps all columns.
  std::vector<bool> activeMask(blocksize, true);

  double smallestResidualNorm = std::numeric_limits<double>::max();
  int64_t iterationNumber = -1;
  bool restart = true;        // no valid P yet
  bool forcedRestart = false; // set when residuals blow up
  bool explicitGramFlag = false;
  std::string criterion = "maxiterations";

  arma::vec residualNorms;

  while (iterationNumber < max_iterations) {
    ++iterationNumber;

    // Residuals R = A X - X Lambda and their column norms.
    R = AX - X * arma::diagmat(lambda);
    residualNorms = column_norms(R, dot);
    double residualNorm = arma::sum(residualNorms) / blocksize;

    lambdaHistory.push_back(lambda.t());
    residualHistory.push_back(residualNorms.t());

    // Track the best (smallest mean-residual) iterate to return.
    if (residualNorm < smallestResidualNorm) {
      smallestResidualNorm = residualNorm;
      bestX = X;
      bestLambda = lambda;
    } else if (residualNorm >
               math::ipow(2, restart_control) * smallestResidualNorm) {
      // Divergence guard: recompute A X exactly and force a P restart.
      forcedRestart = true;
      multiplyA(X, AX);
    }

    Log(1, "LOBPCG iteration {}: mean residual {:.2e}", iterationNumber,
        residualNorm);

    // Converged once the lowest neigs columns are all below tolerance.
    if (residualNorms.head(neigs).max() < tol) {
      criterion = "converged";
      break;
    }

    // Update the active (unconverged) set, monotonically.
    std::vector<arma::uword> active;
    for (int64_t i = 0; i < blocksize; ++i) {
      if (activeMask[i] && (residualNorms(i) <= tol)) {
        activeMask[i] = false;
      }
      if (activeMask[i]) {
        active.push_back((arma::uword)i);
      }
    }
    int64_t currentBlockSize = (int64_t)active.size();
    Log(2, "  active (unconverged) block size {} / {}", currentBlockSize,
        blocksize);
    if (currentBlockSize == 0) {
      criterion = "converged";
      break;
    }
    arma::uvec activeIdx(active);

    // Active residuals, orthogonalized against X (keeps the trial subspace
    // [X, R, P] from re-exploring the X subspace, improving conditioning),
    // then orthonormalized.
    mat_t aR = R.cols(activeIdx);
    aR -= X * dot(X, aR);
    if (!orthonormalize(aR, invR, dot)) {
      Log.warn("LOBPCG: residual orthonormalization failed at iteration {} "
               "(mean residual {:.2e}, tol {:.2e}); returning best iterate.",
               iterationNumber, residualNorm, tol);
      criterion = "failed";
      break;
    }
    mat_t aAR(n, currentBlockSize);
    multiplyA(aR, aAR); // the single MVM per iteration

    // Prepare the conjugate directions P (implicit form, see below).
    mat_t aP, aAP;
    if (iterationNumber > 0) {
      aP = P.cols(activeIdx);
      aAP = AP.cols(activeIdx);
      mat_t invP;
      if (orthonormalize(aP, invP, dot)) {
        aAP = aAP * invP; // keep A P consistent with the normalization
        restart = forcedRestart;
      } else {
        restart = true; // P degenerate -> drop it this iteration
      }
    }

    // Near convergence the assumption of exact orthonormality is swamped by
    // roundoff, so switch to explicitly formed, symmetrized Gram matrices.
    // Monotone latch: once explicit, always explicit.
    double myeps = std::sqrt(std::numeric_limits<double>::epsilon());
    if (residualNorms.max() <= myeps) {
      explicitGramFlag = true;
    }

    // Common blocks (A-weighted).
    mat_t gramXAR = dot(X, aAR);
    mat_t gramRAR = dot(aR, aAR);

    // Overlap blocks; exact-orthonormality shortcuts when far from convergence.
    mat_t gramXX, gramRR, gramXR;
    if (explicitGramFlag) {
      gramRAR = (gramRAR + gramRAR.t()) / 2.0;
      gramXAX = dot(X, AX);
      gramXAX = (gramXAX + gramXAX.t()) / 2.0;
      gramXX = dot(X, X);
      gramRR = dot(aR, aR);
      gramXR = dot(X, aR);
    } else {
      gramXAX = arma::diagmat(arma::conv_to<arma::Col<coeff_t>>::from(lambda));
      gramXX = arma::eye<mat_t>(blocksize, blocksize);
      gramRR = arma::eye<mat_t>(currentBlockSize, currentBlockSize);
      gramXR = mat_t(blocksize, currentBlockSize, arma::fill::zeros);
    }

    mat_t gramA, gramB;
    bool have_solution = false;
    if (!restart) {
      mat_t gramXAP = dot(X, aAP);
      mat_t gramRAP = dot(aR, aAP);
      mat_t gramPAP = dot(aP, aAP);
      mat_t gramXP = dot(X, aP);
      mat_t gramRP = dot(aR, aP);
      mat_t gramPP;
      if (explicitGramFlag) {
        gramPAP = (gramPAP + gramPAP.t()) / 2.0;
        gramPP = dot(aP, aP);
      } else {
        gramPP = arma::eye<mat_t>(currentBlockSize, currentBlockSize);
      }
      gramA =
          sym_block_3x3(gramXAX, gramXAR, gramXAP, gramRAR, gramRAP, gramPAP);
      gramB = sym_block_3x3(gramXX, gramXR, gramXP, gramRR, gramRP, gramPP);
      have_solution = eig_sym_gen(lambda, coeff, gramA, gramB);
      if (!have_solution) {
        Log.warn("LOBPCG: generalized eigensolve failed at iteration {}, "
                 "restarting without conjugate directions.",
                 iterationNumber);
        restart = true; // fall back to the [X, R] subspace below
      }
    }

    if (!have_solution) { // restart path: trial subspace [X, R] only
      gramA = sym_block_2x2(gramXAX, gramXAR, gramRAR);
      gramB = sym_block_2x2(gramXX, gramXR, gramRR);
      if (!eig_sym_gen(lambda, coeff, gramA, gramB)) {
        Log.warn("LOBPCG: generalized eigensolve failed at iteration {} on the "
                 "restart subspace; returning best iterate.",
                 iterationNumber);
        criterion = "failed";
        break;
      }
    }

    // Keep the blocksize smallest Ritz pairs (eig_sym_gen returns ascending):
    // select the lowest-blocksize eigenvectors (columns), then split their rows
    // into the X, R (, P) components of the trial subspace.
    lambda = arma::vec(lambda.head(blocksize));
    coeff = mat_t(coeff.head_cols(blocksize));
    mat_t cX = coeff.rows(0, blocksize - 1);
    mat_t cR = coeff.rows(blocksize, blocksize + currentBlockSize - 1);

    // New implicit conjugate direction: pp is the part of the new Ritz vector
    // coming from R (and old P), *excluding* X. Storing this combination rather
    // than the difference X_new - X_old avoids catastrophic cancellation near
    // convergence -- the central stabilizing trick of LOBPCG.
    mat_t pp = aR * cR;
    mat_t app = aAR * cR;
    if (!restart) {
      mat_t cP = coeff.rows(blocksize + currentBlockSize,
                            blocksize + 2 * currentBlockSize - 1);
      pp += aP * cP;
      app += aAP * cP;
    }

    X = X * cX + pp;
    AX = AX * cX + app; // updated by linear combination -> no extra MVM
    P = pp;
    AP = app;
  } // main iteration loop

  // Final exact Rayleigh-Ritz on the best iterate: clean up eigenpairs and make
  // the returned vectors orthonormal Ritz vectors of A.
  X = bestX;
  multiplyA(X, AX);
  gramXAX = dot(X, AX);
  gramXAX = (gramXAX + gramXAX.t()) / 2.0;
  mat_t gramXX = dot(X, X);
  gramXX = (gramXX + gramXX.t()) / 2.0;
  if (!eig_sym_gen(lambda, coeff, gramXAX, gramXX)) {
    // Fall back to the pre-postprocessing best result.
    lambda = bestLambda;
  } else {
    X = X * coeff;
    AX = AX * coeff;
  }
  R = AX - X * arma::diagmat(lambda);
  residualNorms = column_norms(R, dot);

  if (residualNorms.head(neigs).max() > tol && criterion == "maxiterations") {
    Log.warn("LOBPCG: reached max_iterations={} with accuracy {:.2e} > tol "
             "{:.2e} for the lowest {} eigenvalue(s).",
             max_iterations, residualNorms.head(neigs).max(), tol, neigs);
  }

  arma::mat lambdaHist(lambdaHistory.size(), blocksize);
  arma::mat residualHist(residualHistory.size(), blocksize);
  for (size_t i = 0; i < lambdaHistory.size(); ++i) {
    lambdaHist.row(i) = lambdaHistory[i];
    residualHist.row(i) = residualHistory[i];
  }

  return {lambda,    residualNorms, iterationNumber + 1,
          criterion, lambdaHist,    residualHist};
}
XDIAG_CATCH

} // namespace xdiag::linalg
