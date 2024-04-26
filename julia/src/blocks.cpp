#include "operators.hpp"

namespace xdiag::julia {

void define_spinhalf(jlcxx::Module &mod) {
  mod.add_type<Spinhalf>("Spinhalf")
      .constructor<int64_t>()
      .constructor<int64_t, int64_t>()
      .method("n_sites", &Spinhalf::n_sites)
      .method("size", &Spinhalf::size)
      .method("isreal", &Spinhalf::isreal)
      .method("iscomplex", &Spinhalf::iscomplex);
}

void define_tj(jlcxx::Module &mod) {
  mod.add_type<tJ>("tJ")
      .constructor<int64_t, int64_t, int64_t>()
      .method("n_sites", &tJ::n_sites)
      .method("n_up", &tJ::n_up)
      .method("n_dn", &tJ::n_dn)
      .method("size", &tJ::size)
      .method("isreal", &tJ::isreal)
      .method("iscomplex", &tJ::iscomplex);
}

void define_electron(jlcxx::Module &mod) {
  mod.add_type<Electron>("Electron")
      .constructor<int64_t>()
      .constructor<int64_t, int64_t, int64_t>()
      .method("n_sites", &Electron::n_sites)
      .method("n_up", &Electron::n_up)
      .method("n_dn", &Electron::n_dn)
      .method("size", &Electron::size)
      .method("isreal", &Electron::isreal)
      .method("iscomplex", &Electron::iscomplex);
}
} // namespace xdiag::julia
