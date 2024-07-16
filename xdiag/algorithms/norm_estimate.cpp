#include "norm_estimate.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/random/hash.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag {

// Abstract implementation

//
// Implementing the algorithm 4.1 described in
// "FORTRAN codes for estimating the
// one-norm of a real or complex matrix, with applications to condition
// estimation"
//
// Nicholas J. Higham
// https://dl.acm.org/doi/abs/10.1145/50063.214386
//

template <typename coeff_t, typename apply_A_f, typename apply_A_T_f,
          typename norm_f, typename norm1_f, typename norminf_f>
double norm_estimate(apply_A_f &&apply_A, apply_A_T_f &&apply_A_T,
                     norm_f &&norm, norm1_f &norm1, norminf_f &norminf,
                     int64_t dim, int64_t size, int64_t n_max_attempts,
                     uint64_t seed) {
  // apply_A: function to apply an operator to a vector
  // apply_A_T: function to apply the transpose of that operator to a vector
  // norm: function to compute the 2-norm
  // norm1: function to compute the 1-norm
  // norminf: function to compute the inf-norm
  // dim: dimension of the full vector
  // size: size of the local vector (only different from dim when distributed)
  // n_max_attempts: number of attempts to compute norm
  // seed: random seed
  using vec_t = arma::Col<coeff_t>;

  arma::arma_rng::set_seed(seed);
  vec_t e = vec_t(size, arma::fill::randn);
  e /= norm(e);
  vec_t v = apply_A(e);

  // dim need not necessarily be size for distributed
  if ((dim == 1) && (size == 1)) {
    return std::abs(v(0));
  }

  double gamma = norm1(v);
  vec_t xsi = arma::sign(v);

  vec_t x = apply_A_T(xsi);

  for (int k = 2; k <= n_max_attempts; ++k) {

    // j = min {i ; |x_i| = || x ||_inf}
    double xnorm = norminf(x);

    // v = A e_j
    e.zeros();
    int64_t j = 0;
    for (; j < size; ++j) {
      if (std::abs(std::abs(x(j)) - xnorm) < 1e-14) {
        e(j) = 1.0;
        break;
      }
    }
    v = apply_A(e);

    double gamma_bar = gamma;
    gamma = norm1(v);

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
  for (int64_t i = 1; i <= size; ++i) {
    if (i % 2 == 0) {
      x(i - 1) = 1.0 * (1.0 + (double)(i - 1) / (double)(size - 1));
    } else {
      x(i - 1) = -1.0 * (1.0 + (double)(i - 1) / (double)(size - 1));
    }
  }
  x = apply_A(x);

  double xnorm = norm1(x);
  if (2 * xnorm / (3 * dim) > gamma) {
    gamma = 2 * xnorm / (3 * dim);
  }

  return gamma;
}

double norm_estimate(OpSum const &ops, block_variant_t const &block,
                     int64_t n_max_attempts, uint64_t seed) try {
  return std::visit(
      [&](auto &&block) {
        return norm_estimate(ops, block, n_max_attempts, seed);
      },
      block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return 0.;
}

template <typename block_t>
double norm_estimate(OpSum const &ops, block_t const &block,
                     int64_t n_max_attempts, uint64_t seed) try {
  int iter = 1;
  auto apply_A = [&iter, &ops, &block](arma::cx_vec const &v) {
    auto ta = rightnow();
    auto w = arma::cx_vec(v.n_rows, arma::fill::zeros);
    apply(ops, block, v, block, w);
    Log(2, "Norm estimation iteration {}", iter);
    timing(ta, rightnow(), "MVM", 2);
    ++iter;
    return w;
  };
  auto norm_f = [&block](arma::cx_vec const &v) { return norm(block, v); };
  auto norm1_f = [&block](arma::cx_vec const &v) { return norm1(block, v); };
  auto norminf_f = [&block](arma::cx_vec const &v) {
    return norminf(block, v);
  };
#ifdef XDIAG_USE_MPI
  if (isdistributed(block)) {
    seed = random::hash_combine(random::hash(block), seed);
  }
#endif

  return norm_estimate<complex>(apply_A, apply_A, norm_f, norm1_f, norminf_f,
                                block.dim(), block.size(), n_max_attempts,
                                seed);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return 0.;
}

template double norm_estimate(OpSum const &ops, Spinhalf const &block, int64_t,
                              uint64_t);
template double norm_estimate(OpSum const &ops, tJ const &block, int64_t,
                              uint64_t);
template double norm_estimate(OpSum const &ops, Electron const &block, int64_t,
                              uint64_t);
#ifdef XDIAG_USE_MPI
template double norm_estimate(OpSum const &ops, tJDistributed const &block,
                              int64_t, uint64_t);

#endif

template <typename coeff_t>
double norm_estimate(arma::Mat<coeff_t> const &A, int64_t n_max_attempts,
                     uint64_t seed) try {
  auto apply_A = [&A](arma::Col<coeff_t> const &v) { return A * v; };
  auto apply_A_T = [&A](arma::Col<coeff_t> const &v) {
    arma::Col<coeff_t> vv = A.t() * v;
    return vv;
  };
  auto norm_f = [](arma::Col<coeff_t> const &v) { return arma::norm(v); };
  auto norm1_f = [](arma::Col<coeff_t> const &v) { return arma::norm(v, 1); };
  auto norminf_f = [](arma::Col<coeff_t> const &v) {
    return arma::norm(v, "inf");
  };
  return norm_estimate<coeff_t>(apply_A, apply_A_T, norm_f, norm1_f, norminf_f,
                                A.n_cols, A.n_cols, n_max_attempts, seed);

} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return 0.;
}

template double norm_estimate(arma::mat const &, int64_t, uint64_t);
template double norm_estimate(arma::cx_mat const &, int64_t, uint64_t);
} // namespace xdiag
