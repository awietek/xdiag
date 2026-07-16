#include "basis.hpp"

namespace xdiag::basis {

std::size_t create_basis_type_id() {
  static std::size_t counter = 0;
  return counter++;
}
  
} // namespace xdiag::basis
