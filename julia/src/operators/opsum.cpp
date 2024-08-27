#include "opsum.hpp"

namespace xdiag::julia {

void define_opsum(jlcxx::Module &mod) {

  mod.add_type<OpSum>("cxx_OpSum")
      .constructor<>()
      .constructor<VectorOp const &>()
      .method("size",
              [](OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(ops.size()) })
      .method("defined",
              [](OpSum const &ops, std::string name) {
                JULIA_XDIAG_CALL_RETURN(ops.defined(name))
              })
      .method("getindex",
              [](OpSum const &ops, std::string name) {
                JULIA_XDIAG_CALL_RETURN(ops[name])
              })
      .method("setindex!",
              [](OpSum &ops, Coupling const &cpl, std::string name) {
                JULIA_XDIAG_CALL_VOID(ops[name] = cpl)
              })
      .method("couplings",
              [](OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(ops.couplings()) })
      .method("isreal",
              [](OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(ops.isreal()) })
      .method(
          "isexplicit",
          [](OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(ops.isexplicit()) })
      .method("+", [](OpSum const &ops,
                      Op const &op) { JULIA_XDIAG_CALL_RETURN(ops + op) })
      .method("+", [](OpSum const &ops, OpSum const &ops2) {
        JULIA_XDIAG_CALL_RETURN(ops + ops2)
      });
  // .method("begin", &OpSum::begin)
  // .method("end", &OpSum::end);

  mod.method("to_string", [](OpSum const &ops) { return to_string(ops); });
}
} // namespace xdiag::julia
