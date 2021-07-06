#include "tmatrix.h"

#include <hydra/utils/complex.h>

namespace hydra {
template <class coeff_t>
void Tmatrix<coeff_t>::append(real_type alpha, real_type beta) {
  alphas_.push_back(alpha);
  betas_.push_back(beta);
}

template <class coeff_t>
lila::Vector<lila::real_t<coeff_t>> Tmatrix<coeff_t>::eigenvalues() const {
  if (size() == 0)
    return lila::Vector<lila::real_t<coeff_t>>();
  else
    return lila::EigenvaluesSymTridiag(alphas_, betas_);
}

template <class coeff_t>
std::pair<lila::Vector<lila::real_t<coeff_t>>,
          lila::Matrix<lila::real_t<coeff_t>>>
Tmatrix<coeff_t>::eigen() const {
  return lila::EigenSymTridiag(alphas_, betas_);
}

template class Tmatrix<double>;
template class Tmatrix<complex>;

} // namespace hydra
