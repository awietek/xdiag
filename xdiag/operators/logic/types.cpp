#include "types.hpp"

#include <algorithm>

namespace xdiag {
bool is_known_type(std::string type) {
  return std::find(known_types.begin(), known_types.end(), type) !=
         known_types.end()
}
bool is_real_type(std::string type) {
  return std::find(real_types.begin(), real_types.end(), type) !=
         real_types.end()
}

bool is_cplx_type(std::string type) {
  return std::find(cplx_types.begin(), cplx_types.end(), type) !=
         cplx_types.end()
}

} // namespace xdiag
