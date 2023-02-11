#include "norm_estimate.h"

namespace hydra {

double norm_estimate(std::function<arma::vec(arma::vec const &)> const &apply_A,
                     idx_t N, int n_max_attempts) {
  using namespace arma;

  vec x(N, arma::fill::ones);
  x /= (double)N;

  vec y(N, arma::fill::zeros);
  vec xsi(N, arma::fill::zeros);
  vec z(N, arma::fill::zeros);

  vec v(N);
  for (idx_t i = 0; i < N; ++i) {
    if (i % 2 == 0) {
      v(i) = 1.0 + (double)i / (double)(N-1);
    } else {
      v(i) = -1.0 - (double)i / (double)(N-1);
    }
  }

  for (int attempt = 0; attempt < n_max_attempts; ++attempt) {
    y = apply_A(x);
    xsi = sign(y);
    z = apply_A(xsi);
    double znorm = norm(z, "inf");
    if (znorm <= dot(z, x) && attempt >= 2) {
      double nrm1 = norm(y, 1);
      auto w = apply_A(v);
      double nrm2 = norm(w, 1) / norm(v, 1);
      double gamma = std::max(nrm1, nrm2);
      Log(2, "Operator norm: {} , estimated in {} steps", gamma, attempt);
      return gamma;
    }
    x.zeros();

    // Compute maximizing index j
    idx_t j = 0;
    for (; j < N; ++j) {
      if (std::abs(std::abs(z(j)) - znorm) < 1e-14) {
        break;
      }
    }
    x(j) = 1.0;
  }
  Log.err("Error in norm_estimate: algorithm exceeded maximum number ({}) of "
          "attempts",
          n_max_attempts);
  return 0.;
}

} // namespace hydra
