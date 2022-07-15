#include "complex.h"

#include <hydra/utils/logger.h>

namespace hydra::utils {

void warn_if_complex(complex x, std::string msg) {
  if (!lila::close(lila::imag(x), 0.)) {
    Log.warn("WarningComplexNumber: {}", msg);
  }
}

} // namespace hydra::utils
