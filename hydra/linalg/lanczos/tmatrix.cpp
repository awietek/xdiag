#include "tmatrix.h"

#include <hydra/utils/complex.h>

namespace hydra {

void Tmatrix::append(double alpha, double beta) {
  alphas_.push_back(alpha);
  betas_.push_back(beta);
}
void Tmatrix::pop() {
  alphas_.pop_back();
  betas_.pop_back();
}

lila::Vector<double> Tmatrix::eigenvalues() const {
  if (size() == 0)
    return lila::Vector<double>();
  else
    return lila::EigenvaluesSymTridiag(alphas_, betas_);
}

std::pair<lila::Vector<double>, lila::Matrix<double>> Tmatrix::eigen() const {
  return lila::EigenSymTridiag(alphas_, betas_);
}

} // namespace hydra
