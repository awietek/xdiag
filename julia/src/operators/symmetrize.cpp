#include "symmetrize.hpp"

namespace xdiag::julia {

void define_symmetrize(jlcxx::Module &mod) {

  mod.method("cxx_symmetrize", [](Op const &op, PermutationGroup const &group) {
    JULIA_XDIAG_CALL_RETURN(symmetrize(op, group));
  });
  mod.method("cxx_symmetrize",
             [](OpSum const &ops, PermutationGroup const &group) {
               JULIA_XDIAG_CALL_RETURN(symmetrize(ops, group));
             });
  mod.method("cxx_symmetrize", [](Op const &op, PermutationGroup const &group,
                                  Representation const &irrep) {
    JULIA_XDIAG_CALL_RETURN(symmetrize(op, group, irrep));
  });

  mod.method("cxx_symmetrize",
             [](OpSum const &ops, PermutationGroup const &group,
                Representation const &irrep) {
               JULIA_XDIAG_CALL_RETURN(symmetrize(ops, group, irrep));
             });
}

} // namespace xdiag::julia
