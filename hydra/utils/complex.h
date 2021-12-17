#pragma once

#include <hydra/common.h>
#include <lila/all.h>

namespace hydra::utils {

template <class coeff_t>
coeff_t ForceReal(lila::real_t<coeff_t> number, bool warn, std::string msg) {
  (void)warn;
  (void)msg;
  return (coeff_t)number;
}

template <class coeff_t>
coeff_t ForceReal(lila::complex_t<coeff_t> number, bool warn, std::string msg) {
  if ((warn) && !lila::close(lila::imag(number), 0.))
    std::cerr << msg << std::endl;
  return lila::real(number);
}

void warn_if_complex(complex x, std::string msg = "");

} // namespace hydra::utils
