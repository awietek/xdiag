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

} // namespace xdiag::julia
