// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "permutation_group.hpp"

namespace xdiag::julia {

void define_permutation_group(jlcxx::Module &mod) {

  mod.add_type<PermutationGroup>("cxx_PermutationGroup")
      .method("nsites",
              [](PermutationGroup const &p) {
                JULIA_XDIAG_CALL_RETURN(p.nsites());
              })
      .method("size", [](PermutationGroup const &p) {
        JULIA_XDIAG_CALL_RETURN(p.size());
      });

  mod.method("==", [](PermutationGroup const &a, PermutationGroup const &b) {
    JULIA_XDIAG_CALL_RETURN(a == b);
  });
  mod.method("to_string", [](PermutationGroup const &p) {
    JULIA_XDIAG_CALL_RETURN(to_string(p));
  });

  mod.method("construct_PermutationGroup",
             []() { JULIA_XDIAG_CALL_RETURN(PermutationGroup()); });
  mod.method("construct_PermutationGroup",
             [](int64_t *ptr, int64_t n_permutations, int64_t nsites) {
               JULIA_XDIAG_CALL_RETURN(
                   PermutationGroup(ptr, n_permutations, nsites));
             });
}

} // namespace xdiag::julia
