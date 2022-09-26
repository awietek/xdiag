#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

template <class coeff_t>
void test_exp_sym_v(int n, double precision, int n_tests) {
  using namespace hydra;

  // Log.set_verbosity(2);
  
  for (int seed = 1; seed <= n_tests; ++seed) {
    arma::arma_rng::set_seed(seed);
    arma::Mat<coeff_t> A(n, n, arma::fill::randn);
    A += A.t();
    REQUIRE(A.is_hermitian());

    arma::vec eigs;
    arma::eig_sym(eigs, A);
    A -= eigs(0) * arma::eye<arma::Mat<coeff_t>>(arma::size(A));
    
    arma::Mat<coeff_t> r(1, 1, arma::fill::randn);
    coeff_t tau = r(0, 0);

    HydraPrint(A);
    std::cout << tau << std::endl;

    arma::Col<coeff_t> v0(n, arma::fill::randn);
    arma::Col<coeff_t> Av0 = arma::expmat(tau * A) * v0;

    auto multiply = [&A](arma::Col<coeff_t> const &v, arma::Col<coeff_t> &w) {
      w = A * v;
    };

    arma::Col<coeff_t> Av0_lanczos = exp_sym_v(multiply, v0, tau, precision);
   
    double e = arma::norm(Av0_lanczos - Av0);
    Log("N: {} e: {} precision: {}", n, e, precision);
    REQUIRE(e < precision);
  }
}

TEST_CASE("exp_sym_v", "[algorithms]") {
  using namespace hydra;
  test_exp_sym_v<double>(5, 1e-8, 2);
  test_exp_sym_v<complex>(5, 1e-8, 2);
}
