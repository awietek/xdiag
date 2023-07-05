#include "gram_schmidt.h"

namespace hydra {
template <typename coeff_t>
void gram_schmidt_inplace(arma::Mat<coeff_t> &M, int64_t max_col,
                          std::string method, int iterations) {}

template void gram_schmidt_inplace(arma::Mat<double> &, int64_t, std::string,
                                   int);
template void gram_schmidt_inplace(arma::Mat<complex> &, int64_t, std::string,
                                   int);

template <typename coeff_t>
arma::Mat<coeff_t> gram_schmidt(arma::Mat<coeff_t> const &M, int64_t max_col,
                                std::string method, int iterations) {
  arma::Mat<coeff_t> Q = M;
  gram_schmidt_inplace(Q, max_col, method, iterations);
  return Q;
}
template arma::Mat<double> gram_schmidt(arma::Mat<double> const &M,
                                        int64_t max_col, std::string method,
                                        int iterations);
template arma::Mat<complex> gram_schmidt(arma::Mat<complex> const &M,
                                         int64_t max_col, std::string method,
                                         int iterations);
} // namespace hydra
