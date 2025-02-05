#include "permutation.hpp"

namespace xdiag::julia {
void define_permutation(jlcxx::Module &mod) {
  mod.add_type<Permutation>("cxx_Permutation")
      .constructor<>()
      .constructor<int64_t>()
      .constructor<std::vector<int64_t> const &>()
      .method("inverse",
              [](Permutation const &p) { JULIA_XDIAG_CALL_RETURN(inverse(p)); })
      .method("size",
              [](Permutation const &p) { JULIA_XDIAG_CALL_RETURN(p.size()); })
      .method("array",
              [](Permutation const &p) { JULIA_XDIAG_CALL_RETURN(p.array()); });

  mod.method("cxx_multiply", [](Permutation const &p1, Permutation const &p2) {
    JULIA_XDIAG_CALL_RETURN(multiply(p1, p2));
  });

  mod.method("cxx_pow", [](Permutation const &p, int64_t power) {
    JULIA_XDIAG_CALL_RETURN(pow(p, power));
  });

  mod.method("==", [](Permutation const &a, Permutation const &b) {
    JULIA_XDIAG_CALL_RETURN(a == b);
  });
  mod.method("to_string", [](Permutation const &p) {
    JULIA_XDIAG_CALL_RETURN(to_string(p));
  });
}

} // namespace xdiag::julia
