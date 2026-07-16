#include "ipow.hpp"

namespace xdiag::math {

int64_t ipow(int64_t base, int64_t exp) {
  int64_t result = 1;

  while (exp > 0) {
    if (exp & 1) {
      result *= base;
    }
    base *= base;
    exp >>= 1;
  }
  return result;
}

} // namespace xdiag::math
