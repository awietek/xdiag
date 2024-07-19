#include "basis_electron.hpp"

namespace xdiag::basis {

int64_t dim(BasisElectron const &basis) {
  return std::visit([&](auto &&b) { return b.dim(); }, basis);
}
int64_t size(BasisElectron const &basis) {
  return std::visit([&](auto &&b) { return b.size(); }, basis);
}

template <typename bit_t> bool has_bit_t(BasisElectron const &basis) try {
  return std::visit(
      [](auto &&b) {
        using basis_t = typename std::decay<decltype(b)>::type;
        return std::is_same<bit_t, typename basis_t::bit_t>::value;
      },
      basis);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template bool has_bit_t<uint32_t>(BasisElectron const &basis);
template bool has_bit_t<uint64_t>(BasisElectron const &basis);

} // namespace xdiag::basis
