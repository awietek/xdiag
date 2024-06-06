#include "operators.hpp"

namespace xdiag::julia {

void define_spinhalf(jlcxx::Module &mod) {
  mod.add_type<Spinhalf>("cxx_Spinhalf")
      .constructor<int64_t>()
      .constructor<int64_t, int64_t>()
      .constructor<int64_t, PermutationGroup, Representation>()
      .constructor<int64_t, int64_t, PermutationGroup, Representation>()
      .method("n_sites", &Spinhalf::n_sites)
      .method("size", &Spinhalf::size)
      .method("isreal", &Spinhalf::isreal);
}

void define_tj(jlcxx::Module &mod) {
  mod.add_type<tJ>("cxx_tJ")
      .constructor<int64_t, int64_t, int64_t>()
      .constructor<int64_t, int64_t, int64_t, PermutationGroup,
                   Representation>()
      .method("n_sites", &tJ::n_sites)
      .method("n_up", &tJ::n_up)
      .method("n_dn", &tJ::n_dn)
      .method("size", &tJ::size)
      .method("isreal", &tJ::isreal);
}

void define_electron(jlcxx::Module &mod) {
  mod.add_type<Electron>("cxx_Electron")
      .constructor<int64_t>()
      .constructor<int64_t, int64_t, int64_t>()
      .constructor<int64_t, PermutationGroup, Representation>()
      .constructor<int64_t, int64_t, int64_t, PermutationGroup,
                   Representation>()
      .method("n_sites", &Electron::n_sites)
      .method("n_up", &Electron::n_up)
      .method("n_dn", &Electron::n_dn)
      .method("size", &Electron::size)
      .method("isreal", &Electron::isreal);
}
} // namespace xdiag::julia
