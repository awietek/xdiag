#include "norm_estimate.h"

#include <hydra/algebra/algebra.h>
#include <hydra/algebra/apply.h>
#include <hydra/utils/timing.h>

namespace hydra {

double norm_estimate(BondList const &bonds, block_variant_t const &block) try {
  return std::visit([&](auto &&block) { return norm_estimate(bonds, block); },
                    block);
} catch (...) {
  HydraRethrow("Error dispatching operator norm estimate");
  return 0.;
}

template <typename block_t>
double norm_estimate(BondList const &bonds, block_t const &block) try {
  int iter = 1;
  auto apply_A = [&iter, &bonds, &block](arma::cx_vec const &v) {
    auto ta = rightnow();
    auto w = arma::cx_vec(v.n_rows, arma::fill::zeros);
    apply(bonds, block, v, block, w);
    Log(2, "Norm estimation iteration {}", iter);
    timing(ta, rightnow(), "MVM", 2);
    ++iter;
    return w;
  };
  return norm_estimate_cplx(apply_A, apply_A, block.size());
} catch (...) {
  HydraRethrow("Error computing operator norm estimate");
  return 0.;
}

template double norm_estimate(BondList const &bonds, Spinhalf const &block);
template double norm_estimate(BondList const &bonds, tJ const &block);
template double norm_estimate(BondList const &bonds, Electron const &block);

//
// Implementing the algorithm 4.1 described in
// "FORTRAN codes for estimating the
// one-norm of a real or complex matrix, with applications to condition
// estimation"
//
// Nicholas J. Higham
// https://dl.acm.org/doi/abs/10.1145/50063.214386
//

template <typename coeff_t>
double norm_estimate(
    std::function<arma::Col<coeff_t>(arma::Col<coeff_t> const &)> const
        &apply_A,
    std::function<arma::Col<coeff_t>(arma::Col<coeff_t> const &)> const
        &apply_A_T,
    int64_t N, int n_max_attempts, uint64_t seed) {
  using namespace arma;

  using vec_t = Col<coeff_t>;

  arma_rng::set_seed(seed);
  vec_t e = vec_t(N, arma::fill::randn);
  e /= norm(e);
  vec_t v = apply_A(e);

  if (N == 1) {
    return std::abs(v(0));
  }

  double gamma = norm(v, 1);
  vec_t xsi = sign(v);
  vec_t x = apply_A_T(xsi);

  for (int k = 2; k <= n_max_attempts; ++k) {

    // j = min {i ; |x_i| = || x ||_inf}
    double xnorm = norm(x, "inf");

    // v = A e_j
    e.zeros();
    int64_t j = 0;
    for (; j < N; ++j) {
      if (std::abs(std::abs(x(j)) - xnorm) < 1e-14) {
        e(j) = 1.0;
        break;
      }
    }
    v = apply_A(e);

    double gamma_bar = gamma;
    gamma = norm(v, 1);
    if (approx_equal(sign(v), xsi, "both", 1e-12, 1e-12) ||
        gamma <= gamma_bar) {
      // Log("k: {}", k);
      break;
    }

    auto xsi = sign(v);
    auto x = apply_A_T(xsi);
  }

  // x = (-1)^ ...
  x.zeros();
  for (int64_t i = 1; i <= N; ++i) {
    if (i % 2 == 0) {
      x(i - 1) = 1.0 * (1.0 + (double)(i - 1) / (double)(N - 1));
    } else {
      x(i - 1) = -1.0 * (1.0 + (double)(i - 1) / (double)(N - 1));
    }
  }
  x = apply_A(x);

  double xnorm = norm(x, 1);
  if (2 * xnorm / (3 * N) > gamma) {
    gamma = 2 * xnorm / (3 * N);
  }

  return gamma;
}

template double norm_estimate(
    std::function<arma::Col<double>(arma::Col<double> const &)> const &apply_A,
    std::function<arma::Col<double>(arma::Col<double> const &)> const
        &apply_A_T,
    int64_t N, int n_max_attempts, uint64_t seed);

template double norm_estimate(
    std::function<arma::Col<complex>(arma::Col<complex> const &)> const
        &apply_A,
    std::function<arma::Col<complex>(arma::Col<complex> const &)> const
        &apply_A_T,
    int64_t N, int n_max_attempts, uint64_t seed);

double
norm_estimate_real(std::function<arma::vec(arma::vec const &)> const &apply_A,
                   std::function<arma::vec(arma::vec const &)> const &apply_A_T,
                   int64_t N, int n_max_attempts, uint64_t seed) {
  return norm_estimate<double>(apply_A, apply_A_T, N, n_max_attempts, seed);
}

double norm_estimate_cplx(
    std::function<arma::cx_vec(arma::cx_vec const &)> const &apply_A,
    std::function<arma::cx_vec(arma::cx_vec const &)> const &apply_A_T,
    int64_t N, int n_max_attempts, uint64_t seed) {
  return norm_estimate<complex>(apply_A, apply_A_T, N, n_max_attempts, seed);
}

} // namespace hydra
