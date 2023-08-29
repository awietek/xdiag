#include "error.h"

#include <iomanip> // std::setw

namespace hydra {
void traceback(const std::exception &e, std::size_t depth) {
  if (depth == 0) {
    std::cerr << "Hydra exception trace:\n";
  }
  std::cerr << "[" << std::setw(3) << depth + 1 << "]: " << e.what() << '\n';
  try {
    std::rethrow_if_nested(e);
  } catch (const std::exception &nested) {
    traceback(nested, depth + 1);
  }
  std::cerr << "\n";
}
} // namespace hydra
