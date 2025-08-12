#pragma once

#include <functional>
#include <vector>
#include <xdiag/extern/armadillo/armadillo>

namespace xdiag::lobpcg {

template <typename coeff_t> class MaskedMatrix {
public:
  MaskedMatrix(arma::Mat<coeff_t> const &m, std::vector<bool> const &mask);
  std::vector<arma::Mat<coeff_t>> const &matrices() const;

  void orthonormalize();
  void orthogonalize_to(arma::Mat<coeff_t> const& X);
private:
  std::vector<arma::Mat<coeff_t>> matrices_;
}

// A * X -> Y
template <typename coeff_t>
void apply(
    std::function<void(arma::Mat<coeff_t> const &, arma::Mat<coeff_t>)> A,
    MaskedMatrix<coeff_t> const &X, MaskedMatrix<coeff_t> &Y);

} // namespace xdiag::lobpcg
