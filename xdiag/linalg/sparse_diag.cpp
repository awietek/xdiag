// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "sparse_diag.hpp"

#include <xdiag/linalg/lanczos/eigs_lanczos.hpp>
#include <xdiag/linalg/lanczos/eigvals_lanczos.hpp>
#include <xdiag/linalg/lobpcg/eigs_lobpcg.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag {

template <typename op_t>
static double eigval0(op_t const &ops, Block const &block, double precision,
                      int64_t max_iterations, int64_t random_seed) try {
  if (dim(block) == 0) {
    Log.warn("Warning: block zero dimensional in eigval0");
    return std::nan("");
  }
  auto res = eigvals_lanczos(ops, block, 1, precision, max_iterations, 1e-7,
                             random_seed);

  if (res.eigenvalues.size() == 0) {
    Log.warn("Warning: Tmatrix zero dimensional in eig0_cplx");
    return std::nan("");
  } else {
    return res.eigenvalues(0);
  }
}
XDIAG_CATCH

double eigval0(OpSum const &ops, Block const &block, double precision,
               int64_t max_iterations, int64_t random_seed) try {
  return eigval0<OpSum>(ops, block, precision, max_iterations, random_seed);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
double eigval0(CSRMatrix<idx_t, coeff_t> const &ops, Block const &block,
               double precision, int64_t max_iterations,
               int64_t random_seed) try {
  return eigval0<CSRMatrix<idx_t, coeff_t>>(ops, block, precision,
                                            max_iterations, random_seed);
}
XDIAG_CATCH

template <typename op_t>
static std::tuple<double, State> eig0(op_t const &ops, Block const &block,
                                      double precision, int64_t max_iterations,
                                      int64_t random_seed) try {
  if (dim(block) == 0) {
    Log.warn("Warning: block zero dimensional in eigval0");
    return {std::nan(""), State()};
  }
  auto res =
      eigs_lanczos(ops, block, 1, precision, max_iterations, 1e-7, random_seed);

  if (res.eigenvalues.size() == 0) {
    Log.warn("Warning: Tmatrix zero dimensional in eig0_cplx");
    return {std::nan(""), State()};
  } else {
    return {res.eigenvalues(0), res.eigenvectors};
  }
}
XDIAG_CATCH

std::tuple<double, State> eig0(OpSum const &ops, Block const &block,
                               double precision, int64_t max_iterations,
                               int64_t random_seed) try {
  return eig0<OpSum>(ops, block, precision, max_iterations, random_seed);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
std::tuple<double, State>
eig0(CSRMatrix<idx_t, coeff_t> const &ops, Block const &block, double precision,
     int64_t max_iterations, int64_t random_seed) try {
  return eig0<CSRMatrix<idx_t, coeff_t>>(ops, block, precision, max_iterations,
                                         random_seed);
}
XDIAG_CATCH

// A few oversampling ("guard") vectors: they buffer the requested band from the
// rest of the spectrum, improving convergence and letting a degenerate
// multiplet at the neigs-th eigenvalue be captured with the correct
// multiplicity. For finer control, call eigs_lobpcg directly.
static constexpr int64_t lobpcg_guard = 3;

template <typename op_t>
static arma::vec eigvals(op_t const &ops, Block const &block, int64_t neigs,
                         double precision, int64_t max_iterations,
                         int64_t random_seed) try {
  if (dim(block) == 0) {
    Log.warn("Warning: block zero dimensional in eigvals");
    return arma::vec();
  }
  auto res = eigs_lobpcg(ops, block, neigs, lobpcg_guard,
			 sqrt(precision),  // reduce precision -> too strict
                         max_iterations, random_seed);
  return res.eigenvalues;
}
XDIAG_CATCH

arma::vec eigvals(OpSum const &ops, Block const &block, int64_t neigs,
                  double precision, int64_t max_iterations,
                  int64_t random_seed) try {
  return eigvals<OpSum>(ops, block, neigs, precision, max_iterations,
                        random_seed);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
arma::vec eigvals(CSRMatrix<idx_t, coeff_t> const &ops, Block const &block,
                  int64_t neigs, double precision, int64_t max_iterations,
                  int64_t random_seed) try {
  return eigvals<CSRMatrix<idx_t, coeff_t>>(ops, block, neigs, precision,
                                            max_iterations, random_seed);
}
XDIAG_CATCH

template <typename op_t>
static std::tuple<arma::vec, State>
eigs(op_t const &ops, Block const &block, int64_t neigs, double precision,
     int64_t max_iterations, int64_t random_seed) try {
  if (dim(block) == 0) {
    Log.warn("Warning: block zero dimensional in eigs");
    return {arma::vec(), State()};
  }
  auto res = eigs_lobpcg(ops, block, neigs, lobpcg_guard,
			 sqrt(precision), // reduce precision -> too strict
                         max_iterations, random_seed);
  return {res.eigenvalues, res.eigenvectors};
}
XDIAG_CATCH

std::tuple<arma::vec, State> eigs(OpSum const &ops, Block const &block,
                                  int64_t neigs, double precision,
                                  int64_t max_iterations,
                                  int64_t random_seed) try {
  return eigs<OpSum>(ops, block, neigs, precision, max_iterations, random_seed);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
std::tuple<arma::vec, State>
eigs(CSRMatrix<idx_t, coeff_t> const &ops, Block const &block, int64_t neigs,
     double precision, int64_t max_iterations, int64_t random_seed) try {
  return eigs<CSRMatrix<idx_t, coeff_t>>(ops, block, neigs, precision,
                                         max_iterations, random_seed);
}
XDIAG_CATCH

// Explicit instantiations of all sparse-matrix (idx_t, coeff_t) overloads.
#define XDIAG_INST(IDX, COEFF)                                                 \
  template double eigval0(CSRMatrix<IDX, COEFF> const &, Block const &, double, \
                          int64_t, int64_t);                                   \
  template std::tuple<double, State> eig0(CSRMatrix<IDX, COEFF> const &,       \
                                          Block const &, double, int64_t,      \
                                          int64_t);                            \
  template arma::vec eigvals(CSRMatrix<IDX, COEFF> const &, Block const &,     \
                             int64_t, double, int64_t, int64_t);               \
  template std::tuple<arma::vec, State> eigs(CSRMatrix<IDX, COEFF> const &,    \
                                             Block const &, int64_t, double,   \
                                             int64_t, int64_t);
XDIAG_INST(int32_t, double)
XDIAG_INST(int32_t, complex)
XDIAG_INST(int64_t, double)
XDIAG_INST(int64_t, complex)
#undef XDIAG_INST

} // namespace xdiag
