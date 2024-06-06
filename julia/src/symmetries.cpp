#include "symmetries.hpp"

namespace xdiag::julia {
void define_permutation(jlcxx::Module &mod) {
  mod.add_type<Permutation>("cxx_Permutation")
      .constructor<>()
      .constructor<std::vector<int64_t> const &>()
      .method("inverse", &Permutation::inverse)
      .method("array", &Permutation::array);

  mod.method("multiply", [](Permutation const &p1, Permutation const &p2) {
    JULIA_XDIAG_CALL_RETURN(multiply(p1, p2));
  });
}

void define_permutation_group(jlcxx::Module &mod) {

  mod.add_type<PermutationGroup>("cxx_PermutationGroup")
      .constructor<>()
      .constructor<std::vector<Permutation> const &>()
      .method("inverse", &PermutationGroup::inverse)
      .method("n_sites", &PermutationGroup::n_sites)
      .method("size", &PermutationGroup::size);

}

void define_representation(jlcxx::Module &mod) {

  mod.add_type<Representation>("cxx_Representation")
      .constructor<>()
      .constructor<std::vector<complex> const &>()
      .method("isreal", &Representation::isreal)
      .method("size", &Representation::size);

  mod.method("multiply",
             [](Representation const &r1, Representation const &r2) {
               JULIA_XDIAG_CALL_RETURN(multiply(r1, r2));
             });
}

} // namespace xdiag::julia
