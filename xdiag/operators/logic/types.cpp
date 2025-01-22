#include "types.hpp"

#include <algorithm>
#include <xdiag/common.hpp>

namespace xdiag {
bool is_known_type(std::string type) {
  return std::find(known_types.begin(), known_types.end(), type) !=
         known_types.end();
}
bool is_real_type(std::string type) {
  return std::find(real_types.begin(), real_types.end(), type) !=
         real_types.end();
}

bool is_cplx_type(std::string type) {
  return std::find(cplx_types.begin(), cplx_types.end(), type) !=
         cplx_types.end();
}

int64_t nsites_of_type(std::string type) try {
  auto it = _nsites_of_type.find(type);
  if (it != _nsites_of_type.end()) {
    return it->second;
  } else {
    XDIAG_THROW(fmt::format("Cannot determine number of sites of Op of type "
                            "\"{}\". This type is unknown.",
                            type));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::string known_types_string() {
  std::string str;
  for (auto type : known_types) {
    str += fmt::format("\"{}\", ", type);
  }
  return str.substr(0, str.size() - 2);
}

} // namespace xdiag
