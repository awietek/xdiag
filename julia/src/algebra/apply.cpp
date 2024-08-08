#include "apply.hpp"

namespace xdiag::julia {

void define_apply(jlcxx::Module &mod) {

  mod.method("cxx_apply",
             [](Op const &op, State const &v, State &w, double precision) {
               JULIA_XDIAG_CALL_VOID(apply(op, v, w, precision));
             });

  mod.method("cxx_apply",
             [](OpSum const &ops, State const &v, State &w, double precision) {
               JULIA_XDIAG_CALL_VOID(apply(ops, v, w, precision));
             });
}

} // namespace xdiag::julia
