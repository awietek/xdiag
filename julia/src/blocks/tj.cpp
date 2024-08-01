#include "tj.hpp"

namespace xdiag::julia {

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

  mod.method("to_string", [](tJ const &r) { return to_string(r); });
}

} // namespace xdiag::julia
