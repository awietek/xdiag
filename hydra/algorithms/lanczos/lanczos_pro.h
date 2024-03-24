#pragma once
#include <limits>
#include <random>

#include <hydra/extern/armadillo/armadillo>

#include <hydra/algorithms/lanczos/lanczos_step.h>
#include <hydra/algorithms/lanczos/tmatrix.h>
#include <hydra/common.h>
#include <hydra/utils/logger.h>
#include <hydra/utils/timing.h>

namespace hydra {

template <typename coeff_t> struct lanczos_pro_result {
  arma::vec alphas;
  arma::vec betas;
  arma::mat tmat;
  arma::vec eigenvalues;
  arma::Mat<coeff_t> V;
  int num_iterations;
  std::string criterion;
  int num_reorthogonalizations;
  arma::vec orthogonality_levels;
};

arma::mat compute_omega(arma::vec const &alpha, arma::vec const &beta,
                        double new_beta, int64_t dim, int last_reortho) {
  // Computes the estimate of the orthogonality omega
  // Following Horst D. Simon, The Lanczos algorithm with partial reortho
  // Mathematics of Computation, Jan 1984, 42, 165, pp. 115-142
  // could be shortened to only compute last row of omega

  // the choice of random numbers theta and psi is a bit black magic,
  // here we are following exactly Horst Simons suggestions which turn
  // out to give satisfactory results

  using namespace arma;
  double eps = std::numeric_limits<double>::epsilon();

  std::mt19937_64 generator(42);
  std::normal_distribution<double> dist1(0.0, 0.3);
  std::normal_distribution<double> dist2(0.0, 0.6);
  std::normal_distribution<double> dist3(0.0, 1.5);

  assert(alpha.size() == beta.size());
  int N = alpha.size();
  mat omega = zeros(N + 1, N + 1);
  for (int k = 0; k < N + 1; ++k) {
    omega(k, k) = 1.0;
  }
  for (int k = 1; k < N + 1; ++k) {
    omega(k, k - 1) =
        eps * sqrt((double)dim) * (beta(1) / new_beta) * dist2(generator);
  }

  for (int j = 0; j < N; ++j) {

    if ((j == last_reortho) || (j == last_reortho - 1)) {
      for (int k = 0; k < j; ++k) {
        omega(j + 1, k) = eps * dist3(generator);
      }
    } else {

      for (int k = 0; k < j; ++k) {
        double theta = eps * dist1(generator);

        if (j < N - 1) {
          theta *= (beta(k + 1) + beta(j + 1));
        } else {
          theta *= (beta(k + 1) + new_beta);
        }

        omega(j + 1, k) = beta(k + 1) * omega(j, k + 1) +
                          (alpha(k) - alpha(j)) * omega(j, k) + theta;
        if (k > 0) {
          omega(j + 1, k) += beta(k) * omega(j, k - 1);
        }
        if (j > 0) {
          omega(j + 1, k) -= beta(j) * omega(j - 1, k);
        }

        omega(j + 1, k) /= (j == N - 1) ? new_beta : beta(j + 1);
      }
    }
  }
  return omega;
}

// Generic Lanczos implementation building multiple vectors
template <class coeff_t, class multiply_f, class convergence_f>
lanczos_pro_result<coeff_t>
lanczos_pro(multiply_f mult, arma::Col<coeff_t> &v0, convergence_f converged,
            int max_iterations = 300, double orthogonality_level = 1e-8,
            double deflation_tol = 1e-7,
            bool compute_orthogonality_level = false) {

  using namespace arma;

#ifdef HYDRA_ENABLE_MPI
  auto dot = [](Col<coeff_t> const &v, Col<coeff_t> const &w) -> coeff_t {
    return DotMPI(v, w);
  };
#else
  auto dot = [](Col<coeff_t> const &v, Col<coeff_t> const &w) -> coeff_t {
    return cdot(v, w);
  };
#endif

  auto norm = [&dot](Col<coeff_t> const &v) {
    return std::sqrt(hydra::real(dot(v, v)));
  };

  lanczos_pro_result<coeff_t> res;

  // Zero dimensional problem -> return defaults
  if (v0.size() == 0) {
    res.num_iterations = 0;
    res.criterion = "zerodimensional";
    res.num_reorthogonalizations = 0;
    return res;
  }

  auto tmatrix = Tmatrix();

  // Initialize Lanczos vectors and tmatrix
  auto v1 = v0;
  Col<coeff_t> w(v0.size(), fill::zeros);
  v0.zeros();
  double alpha = 0.;
  double beta = 0.;

  res.num_reorthogonalizations = 0;
  std::vector<double> orthogonality_levels;

  // Normalize start vector or return if norm is zero
  coeff_t v1_norm = norm(v1);
  if (std::abs(v1_norm) > 1e-12) {
    v1 /= v1_norm;
  } else {
    res.criterion = "v0zero";
    return res;
  }

  Mat<coeff_t> &V = res.V;
  try {
    V = Mat<coeff_t>(v0.size(), max_iterations, fill::zeros);
  } catch (...) {
    Log.err("Error in lanczos_pro: Unable to allocate Lanczos vector matrix, "
            "dim=({},{})",
            v0.size(), max_iterations);
  }
  V.col(0) = v1;

  // Main Lanczos loop
  int iteration = 0;
  int last_reortho = 1;
  while (!converged(tmatrix)) {
    Log(1, "Lanczos (PRO) iteration #{}", iteration + 1);

    auto tls0 = rightnow();

    auto tlsi0 = rightnow();

    lanczos_step(v0, v1, w, alpha, beta, mult, dot);
    auto tlsi1 = rightnow();
    timing(tlsi0, tlsi1, "    time MVM", 2);

    if (iteration > 1) {
      Log(2, "    Estimating ortho level ...");
      auto tom0 = rightnow();
      mat omega = compute_omega(tmatrix.alphas(), tmatrix.betas(), beta,
                                v0.size(), last_reortho);
      auto tom1 = rightnow();
      timing(tom0, tom1, "    time ortho estimate", 2);

      double ortho_level_estimate =
          max(abs(omega(iteration, span(0, iteration - 1))));

      if (compute_orthogonality_level) { // should be used only for testing
        Mat<coeff_t> ovlps2 = (V.t() * V);
        Mat<coeff_t> ovlps =
            ovlps2(span(0, iteration - 1), span(0, iteration - 1));
        double ortho_level_computed =
            max(abs(ovlps(iteration - 1, span(0, iteration - 2))));
        Log(2, "    ortho level estimated: {}, computed: {}",
            ortho_level_estimate, ortho_level_computed);
        orthogonality_levels.push_back(ortho_level_computed);
      } else {
        Log(2, "    ortho level estimated: {}", ortho_level_estimate);
        orthogonality_levels.push_back(ortho_level_estimate);
      }

      if (ortho_level_estimate > orthogonality_level) {
        Log(1, "  Performing reorthogonalization");
        auto to0 = rightnow();

        orthogonalize_inplace(v1, V, iteration - 1);
        orthogonalize_inplace(v0, V, iteration - 2);
        beta = norm(v1);
        last_reortho = iteration;
        ++res.num_reorthogonalizations;

        auto to1 = rightnow();
        timing(to0, to1, "  time reortho", 1);
      }
    }
    tmatrix.append(alpha, beta);
    tmatrix.print_log();

    ++iteration;

    // Finish if Lanczos sequence is exhausted
    if (abs(beta) > deflation_tol) {
      v1 /= beta;
    } else {
      res.criterion = "deflation";
      break;
    }

    if (iteration >= max_iterations) {
      res.criterion = "maxiterations";
      break;
    }

    V.col(iteration) = v1;

    auto tls1 = rightnow();
    timing(tls0, tls1, "time Lanczos step", 1);
    Log(2, "");
  }

  if (converged(tmatrix)) {
    res.criterion = "converged";
  }

  res.alphas = tmatrix.alphas();
  res.betas = tmatrix.betas();
  if (iteration == 1) {
    res.tmat = diagmat(res.alphas);
  } else {
    res.tmat = diagmat(res.alphas) +
               diagmat(res.betas.subvec(0, iteration - 2), -1) +
               diagmat(res.betas.subvec(0, iteration - 2), 1);
  }
  if (iteration != max_iterations) {
    res.V = V.head_cols(iteration);
  }
  res.eigenvalues = tmatrix.eigenvalues();
  res.num_iterations = iteration;
  res.orthogonality_levels = arma::vec(orthogonality_levels);

  return res;
}

} // namespace hydra
