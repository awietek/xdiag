// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "lanczos_convergence.hpp"

#include <xdiag/utils/logger.hpp>

namespace xdiag::lanczos {

bool converged_eigenvalues(Tmatrix const &tmat, int n_eigenvalue,
                           double precision) try {
  if (n_eigenvalue < 1) {
    XDIAG_THROW("Number of eigenvalues to converge (parameter: n_eigenvalue) "
                "must be larger or equal to 1");
  }

  int size = tmat.size();
  if (size <= n_eigenvalue + 1)
    return false;
  else {
    if (std::abs(tmat.betas()(size - 1)) < 1e-8)
      return true;

    auto eigs = tmat.eigenvalues();
    auto tmat_previous = tmat;
    tmat_previous.pop();
    auto eigs_previous = tmat_previous.eigenvalues();

    double residue =
        std::abs(eigs(n_eigenvalue - 1) - eigs_previous(n_eigenvalue - 1)) /
        std::abs(eigs(n_eigenvalue - 1));

    Log(2,
        "convergence -> n_eigenvalue: {}, residue: {:.2e}, precision: {:.2e}",
        n_eigenvalue, residue, precision);

    return (residue < precision);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return false;
}

bool converged_time_evolution(Tmatrix const &tmat, complex tau,
                              double precision, double nrm) try {
  int size = tmat.size();
  if (size < 2) {
    return false;
  } else {
    // Lanczos sequence exhausted
    double beta = tmat.betas()(size - 1);
    if (std::abs(beta) < 1e-8) {
      return true;
    }

    // Prepare extended T-matrix for exponentiation
    auto tmatr = tmat.mat();
    auto tmat_ext = arma::mat(size + 2, size + 2);
    for (arma::uword i = 0; i < tmatr.n_rows; ++i) {
      for (arma::uword j = 0; j < tmatr.n_cols; ++j) {
        tmat_ext(i, j) = tmatr(i, j);
      }
    }
    int ext_size = size + 2;
    tmat_ext.resize(ext_size, ext_size);

    tmat_ext(size - 1, size) = beta;
    tmat_ext(size, size - 1) = beta;

    tmat_ext(size - 1, size) = 0.;
    tmat_ext(size + 1, size) = 1.;

    // Exponentiate extended T-matrix
    arma::cx_mat tmat_ext_exp = arma::expmat(tau * tmat_ext);
    double phi1 = std::abs(nrm * tmat_ext_exp(size, 0));
    double phi2 = std::abs(nrm * tmat_ext_exp(size + 1, 0));

    double error;
    if (phi1 > 10 * phi2) {
      error = phi2;
    } else if (phi1 > phi2) {
      error = (phi1 * phi2) / (phi1 - phi2);
    } else {
      error = phi1;
    }
    return (error < precision);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return false;
}

// template <class coeff_t>
// bool ConvergedRitz(const Tmatrix<coeff_t> &tmat, int n_eigenvalue,
//                    coeff_t precision) {
//   int size = tmat.size();
//   if (size < 2)
//     return false;
//   else {
//     // Lanczos sequence exhausted
//     if (close(beta, (coeff_t)0.))
//       return true;
//     if (size <= n_eigenvalue)
//       return false;

//     // Compute all Ritz residuals
//     auto evecs = Eigen(tmat).eigenvectors;
//     bool conv = true;
//     for (int n = 0; n <= n_eigenvalue; ++n) {
//       int idx = size - 1 - n_eigenvalue;
//       double residue = std::abs(evecs(idx, size - 1) * beta);
//       conv &= (residue < precision);
//     }
//     return conv;
//   }
// }

// template <class coeff_t>
// bool ConvergedFixed(const Tmatrix<coeff_t> &tmat, int n_iterations) {
//   int size = tmat.size();
//   return size >= n_iterations;
// }

} // namespace xdiag::lanczos
