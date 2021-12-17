#include "complex.h"

namespace hydra::utils {

void warn_if_complex(complex x, std::string msg) {
  if (!lila::close(lila::imag(x), 0.)) {
    lila::Log.warn("WarningComplexNumber: {}", msg);
  }
}

} // namespace hydra::utils
