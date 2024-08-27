#include "permutation_group.hpp"

namespace xdiag::julia {

void define_permutation_group(jlcxx::Module &mod) {

  mod.add_type<PermutationGroup>("cxx_PermutationGroup")
      .constructor<>()
      .constructor<VectorPermutation const &>()
      .method("inverse", &PermutationGroup::inverse)
      .method("n_sites", &PermutationGroup::n_sites)
      .method("size", &PermutationGroup::size);

  mod.method("to_string",
             [](PermutationGroup const &p) { return to_string(p); });

}

} // namespace xdiag::julia
