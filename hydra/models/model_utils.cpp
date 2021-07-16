#include "model_utils.h"

namespace hydra::utils {

bool coupling_is_zero(Bond const &bond, Couplings const &couplings) {
  std::string coupling = bond.coupling();
  return (!couplings.defined(coupling)) ||
         lila::close(couplings[coupling], (complex)0.);
}
  
} // namespace hydra
