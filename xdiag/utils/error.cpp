#include "error.h"

#include <iomanip> // std::setw
#include <variant>

#ifdef __APPLE__
/// DIRTY HACK to make std::visit work on old MacOS versions
const char* std::bad_variant_access::what() const noexcept {
    return "bad_variant_access";
}
#endif

namespace xdiag {
void traceback(const std::exception &e, std::size_t depth) {
  if (depth == 0) {
    std::cerr << "XDiag exception trace:\n";
  }
  std::cerr << "[" << std::setw(3) << depth + 1 << "]: " << e.what() << '\n';
  try {
    std::rethrow_if_nested(e);
  } catch (const std::exception &nested) {
    traceback(nested, depth + 1);
  }
  if (depth == 0) {
    std::cerr << "\n";
  }
}
} // namespace xdiag
