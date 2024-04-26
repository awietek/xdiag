#include "operators.hpp"

namespace xdiag::julia {

void define_bond(jlcxx::Module &mod) {
  mod.add_type<Bond>("BondCxx")
      .constructor<>()
      .constructor<std::string, std::string, std::vector<int64_t> const &>()
      .constructor<std::string, complex, std::vector<int64_t> const &>()
      .constructor<std::string, double, std::vector<int64_t> const &>()
      .method("isreal", &Bond::isreal)
      .method("iscomplex", &Bond::iscomplex);
}

void define_bondlist(jlcxx::Module &mod) {
  mod.add_type<BondList>("BondListCxx")
      .constructor<>()
      .constructor<std::vector<Bond> const &>()
      .method("set_coupling", &BondList::set_coupling)
      .method("isreal", &BondList::isreal)
      .method("iscomplex", &BondList::iscomplex);
}
} // namespace xdiag::julia
