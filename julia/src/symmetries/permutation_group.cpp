#include "permutation_group.hpp"

namespace xdiag::julia {

void define_permutation_group(jlcxx::Module &mod) {

  mod.add_type<PermutationGroup>("cxx_PermutationGroup")
      .constructor<>()
      .constructor<arma::Mat<int64_t> const &>()
      .method("nsites",
              [](PermutationGroup const &p) {
                JULIA_XDIAG_CALL_RETURN(p.nsites());
              })
      .method("size", [](PermutationGroup const &p) {
        JULIA_XDIAG_CALL_RETURN(p.size());
      });

  mod.method("to_string", [](PermutationGroup const &p) {
    JULIA_XDIAG_CALL_RETURN(to_string(p));
  });
}

} // namespace xdiag::julia
