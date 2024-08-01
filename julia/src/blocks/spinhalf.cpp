#include "spinhalf.hpp"

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

    mod.method("to_string", [](Spinhalf const &r) { return to_string(r); });
}

} // namespace xdiag::julia
