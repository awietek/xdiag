#pragma once

#include <cmath>
#include <xdiag/common.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

namespace xdiag::lobpcg {

template <typename coeff_t>
bool eig_sym_gen(arma::vec &eigvals, arma::Mat<coeff_t> &eigvecs,
                 arma::Mat<coeff_t> const &A, arma::Mat<coeff_t> const &B) try {
  // solve generalized eigenvalue problem gramA * z = lamda * gram * z

  // see https://www.netlib.org/lapack/lug/node54.html how to do it
  // by reducing it to an ordinary eigenvalue problem with cholesky

  double eps = std::sqrt(std::numeric_limits<double>::epsilon());

  if (!A.is_hermitian(eps)) {
    Log.warn("Gram A matrix found not to be hermitian.");
  }
  if (!B.is_hermitian(eps)) {
    Log.warn("Gram matrix found not to be hermitian.");
  }

  if (!B.is_sympd(1e-6)) {
    B.print();
    eig_sym(B).print();
    printf("hello here\n");
    return false;
  }

  // first do a cholesky decomposition of gram
  auto L = chol(B, "lower");
  arma::Mat<coeff_t> Linv = inv(trimatl(L));

  // Compute adjoint matrix C = L^-1 * gramA * L^-T
  arma::Mat<coeff_t> C = Linv * A * Linv.t();
  bool success = eig_sym(eigvals, eigvecs, C);
  if (!success) {
    printf("hello\n");
    return false;
  }

  // retrieve eigenvectors z = L^-T w
  eigvecs = solve(trimatu(L.t()), eigvecs);
  return true;
}
XDIAG_CATCH

template <typename coeff_t>
static arma::Mat<coeff_t>
block_2x2(arma::Mat<coeff_t> const &A, arma::Mat<coeff_t> const &B,
          arma::Mat<coeff_t> const &C, arma::Mat<coeff_t> const &D) try {
  int64_t nr1 = A.n_rows;
  int64_t nr2 = C.n_rows;
  if ((B.n_rows != nr1) || (D.n_rows != nr2)) {
    XDIAG_THROW("Invalid row dimensions");
  }

  int64_t nc1 = A.n_cols;
  int64_t nc2 = B.n_cols;
  if ((C.n_cols != nc1) || (D.n_cols != nc1)) {
    XDIAG_THROW("Invalid column dimensions");
  }

  int64_t nr = nr1 + nr2;
  int64_t nc = nc1 + nc2;

  arma::Mat<coeff_t> M(nr, nc, arma::fill::zeros);
  using arma::span;
  M.submat(span(0, nr1 - 1), span(0, nc1 - 1)) = A;
  M.submat(span(0, nr1 - 1), span(nc1, nc - 1)) = B;

  M.submat(span(nr1, nr - 1), span(0, nc1 - 1)) = C;
  M.submat(span(nr1, nr - 1), span(nc1, nc - 1)) = D;
  return M;
}
XDIAG_CATCH

template <typename coeff_t>
static arma::Mat<coeff_t>
block_3x3(arma::Mat<coeff_t> const &A, arma::Mat<coeff_t> const &B,
          arma::Mat<coeff_t> const &C, arma::Mat<coeff_t> const &D,
          arma::Mat<coeff_t> const &E, arma::Mat<coeff_t> const &F,
          arma::Mat<coeff_t> const &G, arma::Mat<coeff_t> const &H,
          arma::Mat<coeff_t> const &I) try {
  int64_t nr1 = A.n_rows;
  int64_t nr2 = D.n_rows;
  int64_t nr3 = G.n_rows;
  if (((B.n_rows != nr1) || (C.n_rows != nr1)) ||
      ((E.n_rows != nr2) || (F.n_rows != nr2)) ||
      ((H.n_rows != nr3) || (I.n_rows != nr3))) {
    XDIAG_THROW("Invalid row dimensions");
  }

  int64_t nc1 = A.n_cols;
  int64_t nc2 = B.n_cols;
  int64_t nc3 = C.n_cols;
  if (((D.n_cols != nc1) || (G.n_cols != nc1)) ||
      ((E.n_cols != nc2) || (H.n_cols != nc2)) ||
      ((F.n_cols != nc3) || (I.n_cols != nc3))) {
    XDIAG_THROW("Invalid column dimensions");
  }

  int64_t nr = nr1 + nr2 + nr3;
  int64_t nc = nc1 + nc2 + nc3;

  arma::Mat<coeff_t> M(nr, nc, arma::fill::zeros);
  using arma::span;
  M.submat(span(0, nr1 - 1), span(0, nc1 - 1)) = A;
  M.submat(span(0, nr1 - 1), span(nc1, nc1 + nc2 - 1)) = B;
  M.submat(span(0, nr1 - 1), span(nc1 + nc2, nc - 1)) = C;

  M.submat(span(nr1, nr1 + nr2 - 1), span(0, nc1 - 1)) = D;
  M.submat(span(nr1, nr1 + nr2 - 1), span(nc1, nc1 + nc2 - 1)) = E;
  M.submat(span(nr1, nr1 + nr2 - 1), span(nc1 + nc2, nc - 1)) = F;

  M.submat(span(nr1 + nr2, nr - 1), span(0, nc1 - 1)) = G;
  M.submat(span(nr1 + nr2, nr - 1), span(nc1, nc1 + nc2 - 1)) = H;
  M.submat(span(nr1 + nr2, nr - 1), span(nc1 + nc2, nc - 1)) = I;
  return M;
}
XDIAG_CATCH

template <typename coeff_t>
static std::pair<bool, arma::Mat<coeff_t>>
orthonormalize(arma::Mat<coeff_t> &V) try {
  arma::Mat<coeff_t> VBV = V.t() * V;
  if (VBV.is_sympd()) {
    arma::Mat<coeff_t> VBVinv = inv_sympd(VBV);
    V = V * VBVinv;
    return {true, VBVinv};
  } else {
    return {false, arma::Mat<coeff_t>()};
  }
}
XDIAG_CATCH

static std::vector<int64_t> _get_indx(arma::vec const &_lambda, int64_t num,
                                      bool smallest) try {
  if (_lambda.size() < num) {
    XDIAG_THROW("_lambda.size() < num (this is a bug, please report)");
  }

  // get sorting indices of _lambda
  std::vector<int64_t> ii(_lambda.size());
  std::iota(ii.begin(), ii.end(), 0);
  std::stable_sort(ii.begin(), ii.end(), [&](size_t i1, size_t i2) {
    return _lambda(i1) < _lambda(i2);
  });

  if (smallest) {
    return std::vector<int64_t>(ii.begin(), ii.begin() + num);
  } else {
    std::vector<int64_t> last(ii.end() - num, ii.end());
    std::reverse(last.begin(), last.end());
    return last;
  }
}
XDIAG_CATCH

static arma::vec _select(arma::vec const &v,
                         std::vector<int64_t> const &idx) try {
  if (idx.size() > v.size()) {
    XDIAG_THROW("idx.size() > v.size() (this is a bug, please report)");
  }
  arma::vec w(idx.size(), arma::fill::zeros);
  for (int64_t k = 0; k < idx.size(); ++k) {
    w(k) = v(idx[k]);
  }
  return w;
}
XDIAG_CATCH

template <typename coeff_t>
static arma::Mat<coeff_t> _select_columns(arma::Mat<coeff_t> const &M,
                                          std::vector<int64_t> const &idx) try {
  if (idx.size() > M.n_cols) {
    XDIAG_THROW("idx.size() > M.n_cols (this is a bug, please report)");
  }
  arma::Mat<coeff_t> N(M.n_rows, idx.size(), arma::fill::zeros);
  for (int64_t k = 0; k < idx.size(); ++k) {
    N.col(k) = M.col(idx[k]);
  }
  return N;
}
XDIAG_CATCH

template <typename coeff_t, typename multiply_f>
void lobpcg(multiply_f multiplyA, arma::Mat<coeff_t> &X, double tol,
            int64_t maxiter, bool smallest, int64_t restartControl) try {
  using mat_t = arma::Mat<coeff_t>;
  using rvec_t = arma::vec;

  mat_t &blockVectorX = X;
  mat_t bestBlockVectorX;
  try {
    bestBlockVectorX = blockVectorX;
  } catch (...) {
    XDIAG_THROW("Unable to allocate memory for auxiliary block vectors");
  }

  double residualTolerance = tol;
  if (maxiter < 0) {
    maxiter = 20;
  }

  if ((blockVectorX.n_rows == 0) || (blockVectorX.n_cols == 0)) {
    XDIAG_THROW("Input initial matrix X has zero dimension");
  }

  int64_t n = blockVectorX.n_rows;
  int64_t sizeX = blockVectorX.n_cols;

  mat_t lambdaHistory(maxiter + 3, sizeX, arma::fill::zeros);
  mat_t residualNormsHistory(maxiter + 3, sizeX, arma::fill::zeros);

  Log("Solving standard eigenvalue problem\n");
  Log("matrix size {}", n);
  Log("block size  {}\n", sizeX);

  if (n < 5 * sizeX) {
    XDIAG_THROW(fmt::format("The problem size {} is too small for the block "
                            "size {}. Cannot iterate in LOBPCG reliably.",
                            n, sizeX));
  }

  if (residualTolerance <= 0.0) {
    residualTolerance = std::sqrt(std::numeric_limits<double>::epsilon()) * n;
  }

  Log("before init ortho");

  // orthonormalize initial X
  arma::Mat<coeff_t> XX = blockVectorX.t() * blockVectorX;
  if (XX.is_sympd()) {
    arma::Mat<coeff_t> XXinv = inv_sympd(XX);
    blockVectorX = blockVectorX * XXinv;
  } else {
    XDIAG_THROW("Linearly dependent initial approximations");
  }

  Log("before allocation");

  // Allocate memory for all auxiliary block vectors
  mat_t blockVectorAX;
  mat_t blockVectorP;
  mat_t blockVectorAP;
  mat_t blockVectorR;
  mat_t blockVectorAR;
  mat_t aux;
  try {
    blockVectorAX.resize(n, sizeX);
    blockVectorP.resize(n, sizeX);
    blockVectorAP.resize(n, sizeX);
    blockVectorR.resize(n, sizeX);
    blockVectorAR.resize(n, sizeX);
    aux.resize(n, sizeX);
    blockVectorAX.zeros();
    blockVectorP.zeros();
    blockVectorAP.zeros();
    blockVectorR.zeros();
    blockVectorAR.zeros();
    aux.zeros();
  } catch (...) {
    XDIAG_THROW("Unable to allocate memory for auxiliary block vectors");
  }

  // Compute the initial Ritz vectors: solve the eigenproblem.
  // MVM
  multiplyA(blockVectorX, blockVectorAX);

  mat_t gramXAX = blockVectorX.t() * blockVectorAX;
  rvec_t _lambda;
  mat_t eigBlockVector;
  eig_sym(_lambda, eigBlockVector, gramXAX);
  std::vector<int64_t> ii = _get_indx(_lambda, sizeX, smallest);
  _lambda = _select(_lambda, ii);
  lambdaHistory.row(0) = _lambda.t();

  eigBlockVector = _select_columns(eigBlockVector, ii);

  blockVectorX = blockVectorX * eigBlockVector;
  blockVectorAX = blockVectorAX * eigBlockVector;

  // Active index set.
  std::vector<bool> activeMask(sizeX, true);

  double smallestResidualNorm = std::abs(std::numeric_limits<double>::max());
  int64_t iterationNumber = -1;
  int64_t bestIterationNumber = -1;
  bool restart = true;
  bool forcedRestart = false;
  bool explicitGramFlag = true;

  while (iterationNumber < maxiter) {
    iterationNumber += 1;
    Log(1, "iteration {}", iterationNumber);

    aux = blockVectorX * arma::diagmat(_lambda);
    blockVectorR = blockVectorAX - aux;

    // get 2-norm of columns
    // rvec_t residualNorms = arma::vecnorm(blockVectorR, 2, 0);
    rvec_t residualNorms = real(mat_t(blockVectorR.t() * blockVectorR).diag());
    residualNormsHistory.row(iterationNumber) = residualNorms.t();
    double residualNorm = arma::sum(residualNorms) / sizeX;

    if (residualNorm < smallestResidualNorm) {
      smallestResidualNorm = residualNorm;
      bestIterationNumber = iterationNumber;
      bestBlockVectorX = blockVectorX;
    } else if (residualNorm >
               std::pow(2, restartControl) * smallestResidualNorm) {
      forcedRestart = true;

      // MVM
      multiplyA(blockVectorX, blockVectorAX);
    }

    // Update active mask
    for (int64_t i = 0; i < residualNorms.size(); ++i) {
      if (residualNorms[i] < residualTolerance) {
        activeMask[i] = false;
      }
    }
    // int64_t currentBlockSize =
    //     std::accumulate(activeMask.begin(), activeMask.end(), 0);
    int64_t currentBlockSize = sizeX;
    Log(1, "current block size: {}", currentBlockSize);
    if (Log.verbosity() > 0) {
      Log(1, "eigenvalue(s):");
      _lambda.print();
      Log(1, "residual norm(s):");
      residualNorms.print();
      Log(1, "active mask:");
      for (bool m : activeMask) {
        std::cout << (int)m << " ";
      }
      printf("\n");
    }

    if (currentBlockSize == 0) {
      break;
    }

    // Henceforth, only work on eigenvalues which are not converged already
    auto &activeBlockVectorR = blockVectorR;
    auto &activeBlockVectorAR = blockVectorAR;

    // if (iterationNumber > 0) {
    auto &activeBlockVectorP = blockVectorP;
    auto &activeBlockVectorAP = blockVectorAP;
    // }

    // orthogonalize the residuals to X
    activeBlockVectorR =
        activeBlockVectorR -
        (blockVectorX) * (blockVectorX.t() * activeBlockVectorR);

    auto [success, _] = orthonormalize(activeBlockVectorR);
    if (!success) {
      Log.warn(fmt::format("Failed at iteration {} with accuracies {} not "
                           "reaching the requested tolerance {}.",
                           iterationNumber, residualNorm, residualTolerance));
      break;
    }

    // MVM
    multiplyA(activeBlockVectorR, activeBlockVectorAR);

    if (iterationNumber > 0) {
      auto [success, invR] = orthonormalize(activeBlockVectorP);
      if (success) {
        activeBlockVectorAP = activeBlockVectorAP * invR;
        restart = forcedRestart;
      } else {
        restart = true;
      }
    }

    double myeps = std::sqrt(std::numeric_limits<double>::epsilon());
    if ((residualNorms.max() > myeps) && !explicitGramFlag) {
      explicitGramFlag = false;
    } else {
      // Once explicitGramFlag, forever explicitGramFlag.
      explicitGramFlag = true;
    }

    // Common submatrices
    mat_t gramXAX;
    mat_t gramXAR = blockVectorX.t() * activeBlockVectorAR;
    mat_t gramXAP;
    mat_t gramRAR = blockVectorR.t() * activeBlockVectorAR;
    mat_t gramRAP;
    mat_t gramPAP;

    mat_t gramXX;
    mat_t gramXR;
    mat_t gramXP;
    mat_t gramRR;
    mat_t gramRP;
    mat_t gramPP;
    if (explicitGramFlag) {
      gramRAR = (gramRAR + gramRAR.t()) / 2.0;
      gramXAX = blockVectorX.t() * blockVectorAX;
      gramXAX = (gramXAX + gramXAX.t()) / 2.0;
      gramXX = blockVectorX.t() * blockVectorX;
      gramRR = activeBlockVectorR.t() * activeBlockVectorR;
      gramXR = blockVectorX.t() * activeBlockVectorR;
    } else {
      gramXAX = diagmat(_lambda);
      gramXX = arma::eye<mat_t>(sizeX, sizeX);
      gramRR = arma::eye<mat_t>(currentBlockSize, currentBlockSize);
      gramXR = mat_t(currentBlockSize, currentBlockSize, arma::fill::zeros);
    }

    if (!restart) {
      gramXAP = blockVectorX.t() * activeBlockVectorAP;
      gramRAP = activeBlockVectorR.t() * activeBlockVectorAP;
      gramPAP = activeBlockVectorP.t() * activeBlockVectorAP;
      gramXP = blockVectorX.t() * activeBlockVectorP;
      gramRP = activeBlockVectorR.t() * activeBlockVectorP;

      if (explicitGramFlag) {
        gramPAP = (gramPAP + gramPAP.t()) / 2.0;
        gramPP = activeBlockVectorP.t() * activeBlockVectorP;
      } else {
        gramPP = arma::eye<mat_t>(currentBlockSize, currentBlockSize);
      }

      auto gramA =
          block_3x3(gramXAX, gramXAR, gramXAP, mat_t(gramXAR.t()), gramRAR,
                    gramRAP, mat_t(gramXAP.t()), mat_t(gramRAP.t()), gramPAP);
      auto gram =
          block_3x3(gramXX, gramXR, gramXP, mat_t(gramXR.t()), gramRR, gramRP,
                    mat_t(gramXP.t()), mat_t(gramRP.t()), gramPP);

      bool success = eig_sym_gen(_lambda, eigBlockVector, gramA, gram);
      if (!success) {
        Log.warn(1,
                 "eig_sym failed at iteration {} with error causing a restart.",
                 iterationNumber);
        restart = true;
      }

    } else { // if (restart)
      auto gramA = block_2x2(gramXAX, gramXAR, mat_t(gramXAR.t()), gramRAR);
      auto gram = block_2x2(gramXX, gramXR, mat_t(gramXR.t()), gramRR);

      bool success = eig_sym_gen(_lambda, eigBlockVector, gramA, gram);
      if (!success) {
        Log.warn(1, "eig_sym failed at iteration {} with error.",
                 iterationNumber);
        break;
      }
    }

    ii = _get_indx(_lambda, sizeX, smallest);
    _lambda = _select(_lambda, ii);
    eigBlockVector = _select_columns(eigBlockVector, ii);
    lambdaHistory.row(iterationNumber + 1) = _lambda.t();

    auto eigBlockVectorX = eigBlockVector.rows(0, sizeX - 1);
    auto eigBlockVectorR =
        eigBlockVector.rows(sizeX, sizeX + currentBlockSize - 1);

    // Do I really need pp and app, coulnd't that be some other vect
    mat_t pp = activeBlockVectorR * eigBlockVectorR;
    mat_t app = activeBlockVectorAR * eigBlockVectorR;
    if (!restart) {
      auto eigBlockVectorP = eigBlockVector.rows(sizeX + currentBlockSize,
                                                 eigBlockVector.n_rows - 1);
      pp += activeBlockVectorP * eigBlockVectorP;
      app += activeBlockVectorAP * eigBlockVectorP;
    }

    blockVectorX = blockVectorX * eigBlockVectorX + pp;
    blockVectorAX = blockVectorAX * eigBlockVectorX + app;

    blockVectorP = pp;
    blockVectorAP = app;

  } // main iteration loop

  aux = blockVectorX * arma::diagmat(_lambda);
  blockVectorR = blockVectorAX - aux;

  // // get 2-norm of columns
  // residualNorms = arma::vecnorm(blockVectorR, 2, 0);
  // residualNormsHistory.row(iterationNumber + 1) = residualNorms;
  // residualNorm = arma::sum(residualNorms) / sizeX;

  // if (residualNorm < smallestResidualNorm) {
  //   smallestResidualNorm = residualNorm;
  //   bestIterationNumber = iterationNumber + 1;
  //   bestBlockVectorC = blockVectorX;
  // }
}
XDIAG_CATCH

} // namespace xdiag::lobpcg
