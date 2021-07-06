#include "lanczos_convergence.h"

#include <hydra/utils/complex.h>

namespace hydra {

template <class coeff_t>
bool ConvergedEigenvalues(Tmatrix<coeff_t> const &tmat, int n_eigenvalue,
                          lila::real_t<coeff_t> precision) {
  int size = tmat.size();
  auto alphas = tmat.alphas();
  auto betas = tmat.betas();
  if (size < 2)
    return false;
  else {
    auto eigs = tmat.eigenvalues();

    if (lila::close(betas(size - 1), 0.))
      return true;
    if (size <= n_eigenvalue)
      return false;

    auto tmat_previous =
      Tmatrix<coeff_t>(alphas({0, size - 1}), betas({0, size - 1}));
    auto eigs_previous = tmat_previous.eigenvalues();

    double residue =
        std::abs(eigs(n_eigenvalue) - eigs_previous(n_eigenvalue)) /
        std::abs(eigs(n_eigenvalue));
    return (residue < precision);
  }
}

template bool ConvergedEigenvalues<double>(Tmatrix<double> const &, int,
                                           double precision);
template bool ConvergedEigenvalues<complex>(Tmatrix<complex> const &, int,
                                            double precision);

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

// template <class coeff_t, class ccoeff_t>
// bool ConvergedTimeEvolution(const Tmatrix<coeff_t> &tmat, coeff_t beta,
//                             ccoeff_t tau, coeff_t precision, int max_iter,
//                             coeff_t nrm) {
//   int size = tmat.size();
//   if (size < 2)
//     return false;
//   else {
//     // Lanczos sequence exhausted
//     if (close(beta, (coeff_t)0.))
//       return true;

//     // Prepare extended T-matrix for exponentiation
//     auto tmatr = Matrix<coeff_t>(tmat);
//     auto tmat_ext = Zeros<ccoeff_t>(size, size);
//     for (auto i : tmatr.rows())
//       for (auto j : tmatr.cols())
//         tmat_ext(i, j) = (ccoeff_t)tmatr(i, j);

//     int ext_size = size + 2;
//     tmat_ext.resize(ext_size, ext_size);

//     tmat_ext(size - 1, size) = beta;
//     tmat_ext(size, size - 1) = beta;

//     tmat_ext(size - 1, size) = 0.;
//     tmat_ext(size + 1, size) = 1.;

//     // Exponentiate extended T-matrix
//     auto tmat_ext_exp = ExpM(tmat_ext, tau);
//     coeff_t phi1 = std::abs(nrm * tmat_ext_exp(size, 0));
//     coeff_t phi2 = std::abs(nrm * tmat_ext_exp(size + 1, 0));
//     // printf("phi1: %g, phi2: %g\n", phi1, phi2);
//     coeff_t error;
//     if (phi1 > 10 * phi2)
//       error = phi2;
//     else if (phi1 > phi2)
//       error = (phi1 * phi2) / (phi1 - phi2);
//     else
//       error = phi1;
//     if ((error < precision) || (size == max_iter - 1)) {
//       if (size == max_iter - 1)
//         printf("warning: lanczosTevol not converged in %d steps\n",
//         max_iter);
//       return true;
//     } else
//       return false;
//   }
// }

} // namespace hydra
