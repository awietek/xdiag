#include "electron.hpp"

namespace xdiag::julia {

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

  mod.method("to_string", [](Electron const &r) { return to_string(r); });
}
} // namespace xdiag::julia
