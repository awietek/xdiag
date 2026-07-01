// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "eigs_lobpcg.hpp"

#include <algorithm>

#include <xdiag/algebra/ishermitian.hpp>
#include <xdiag/kernels/apply.hpp>
#include <xdiag/kernels/sparse/apply.hpp>
#include <xdiag/linalg/lobpcg/lobpcg.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/math/dot.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/random_state.hpp>

namespace xdiag {

// Builds the matrix-block MVM and (distributed-aware) dot callbacks for a given
// coefficient type, runs the LOBPCG core, and packages the lowest-neigs
// eigenpairs into a State. Shared by the real and complex code paths.
template <typename coeff_t, typename op_t>
static EigsLobpcgResult run_eigs_lobpcg(op_t const &ops, Block const &block,
                                        int64_t neigs, int64_t blocksize,
                                        double tol, int64_t max_iterations,
                                        int64_t random_seed) {

  // Random initial block. Each column is seeded differently so the columns are
  // linearly independent (fill() does not fold the column index into the seed).
  // normalized=false: normalization would compute the norm of the whole
  // multi-column state (unsupported), and LOBPCG orthonormalizes X anyway.
  State X0(block, isreal<coeff_t>(), blocksize);
  for (int64_t col = 0; col < blocksize; ++col) {
    fill(X0, RandomState(random_seed + col, false), col);
  }

  arma::Mat<coeff_t> X;
  if constexpr (isreal<coeff_t>()) {
    X = X0.matrix(false);
  } else {
    X = X0.matrixC(false);
  }

  auto multiplyA = [&ops, &block](arma::Mat<coeff_t> const &V,
                                  arma::Mat<coeff_t> &W) {
    apply(ops, block, V, block, W);
  };
  auto dot = [&block](arma::Mat<coeff_t> const &V,
                      arma::Mat<coeff_t> const &W) {
    return math::matrix_dot(block, V, W);
  };

  linalg::lobpcg_result_t r =
      linalg::lobpcg(multiplyA, dot, X, neigs, tol, max_iterations);

  // Keep the lowest neigs eigenpairs (block is iterated with guard vectors).
  State eigenvectors(block, isreal<coeff_t>(), neigs);
  if constexpr (isreal<coeff_t>()) {
    eigenvectors.matrix(false) = X.head_cols(neigs);
  } else {
    eigenvectors.matrixC(false) = X.head_cols(neigs);
  }

  return {r.eigenvalues.head(neigs),
          r.residual_norms.head(neigs),
          eigenvectors,
          r.niterations,
          r.criterion,
          r.eigenvalue_history,
          r.residual_norms_history};
}

template <typename op_t>
static EigsLobpcgResult
eigs_lobpcg(op_t const &ops, Block const &block, int64_t neigs, int64_t guard,
            double tol, int64_t max_iterations, int64_t random_seed) try {
  int64_t d = dim(block);
  if (d == 0) {
    Log.warn("Warning: block is zero dimensional in eigs_lobpcg");
    return EigsLobpcgResult();
  }
  if (neigs < 1) {
    XDIAG_THROW("Argument \"neigs\" needs to be >= 1");
  }
  if (guard < 0) {
    XDIAG_THROW("Argument \"guard\" needs to be >= 0");
  }
  if (!ishermitian(ops, block)) {
    XDIAG_THROW(
        "Input operator is not hermitian. LOBPCG can only be applied to "
        "hermitian operators.");
  }

  int64_t blocksize = std::min(neigs + guard, d);
  if (blocksize < neigs) {
    XDIAG_THROW("Block size (neigs + guard) exceeds the block dimension");
  }
  if (d < 5 * blocksize) {
    Log.warn("Warning: block dimension {} is small relative to the LOBPCG "
             "block size {} (recommended dim >= 5*blocksize); convergence may "
             "be unreliable -- consider a dense eigensolver.",
             d, blocksize);
  }

  bool real = isreal(ops) && isreal(block);
  if (real) {
    return run_eigs_lobpcg<double>(ops, block, neigs, blocksize, tol,
                                   max_iterations, random_seed);
  } else {
    return run_eigs_lobpcg<complex>(ops, block, neigs, blocksize, tol,
                                    max_iterations, random_seed);
  }
}
XDIAG_CATCH

EigsLobpcgResult eigs_lobpcg(OpSum const &ops, Block const &block,
                             int64_t neigs, int64_t guard, double tol,
                             int64_t max_iterations, int64_t random_seed) try {
  return eigs_lobpcg<OpSum>(ops, block, neigs, guard, tol, max_iterations,
                            random_seed);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
EigsLobpcgResult eigs_lobpcg(CSRMatrix<idx_t, coeff_t> const &ops,
                             Block const &block, int64_t neigs, int64_t guard,
                             double tol, int64_t max_iterations,
                             int64_t random_seed) try {
  return eigs_lobpcg<CSRMatrix<idx_t, coeff_t>>(ops, block, neigs, guard, tol,
                                                max_iterations, random_seed);
}
XDIAG_CATCH

// Template instantiations for every (idx_t, coeff_t) sparse-matrix combination
#define XDIAG_INST(IDX, COEFF)                                                 \
  template EigsLobpcgResult eigs_lobpcg(CSRMatrix<IDX, COEFF> const &,         \
                                        Block const &, int64_t, int64_t,       \
                                        double, int64_t, int64_t);
XDIAG_INST(int32_t, double)
XDIAG_INST(int32_t, complex)
XDIAG_INST(int64_t, double)
XDIAG_INST(int64_t, complex)
#undef XDIAG_INST

} // namespace xdiag
