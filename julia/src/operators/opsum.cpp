#include "opsum.hpp"

namespace xdiag::julia {

void define_opsum(jlcxx::Module &mod) {

  mod.add_type<OpSum>("cxx_OpSum")
      .constructor<>()
      .constructor<Op const &>()
      .constructor<std::string, Op const &>()
      .constructor<double, Op const &>()
      .constructor<complex, Op const &>()
      .method("plain",
              [](OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(ops.plain()); })
      .method("isreal",
              [](OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(isreal(ops)); })
      .method("+", [](OpSum const &ops,
                      Op const &op) { JULIA_XDIAG_CALL_RETURN(ops + op); })
      .method("+",
              [](OpSum const &ops, OpSum const &ops2) {
                JULIA_XDIAG_CALL_RETURN(ops + ops2)
              })
      .method("-", [](OpSum const &ops,
                      Op const &op) { JULIA_XDIAG_CALL_RETURN(ops - op); })
      .method("-",
              [](OpSum const &ops, OpSum const &ops2) {
                JULIA_XDIAG_CALL_RETURN(ops - ops2);
              })
      .method("*", [](OpSum const &ops,
                      double cpl) { JULIA_XDIAG_CALL_RETURN(cpl * ops); })
      .method("*", [](OpSum const &ops,
                      complex cpl) { JULIA_XDIAG_CALL_RETURN(cpl * ops); })
      .method("/", [](OpSum const &ops,
                      double cpl) { JULIA_XDIAG_CALL_RETURN(ops / cpl); })
      .method("/", [](OpSum const &ops,
                      complex cpl) { JULIA_XDIAG_CALL_RETURN(ops / cpl); })
      .method("setindex!",
              [](OpSum &ops, std::string name, double value) {
                JULIA_XDIAG_CALL_VOID(ops[name] = value);
              })
      .method("setindex!",
              [](OpSum &ops, std::string name, complex value) {
                JULIA_XDIAG_CALL_VOID(ops[name] = value);
              })
      .method("constants", [](OpSum const &ops) {
        JULIA_XDIAG_CALL_RETURN(ops.constants());
      });

  mod.method("to_string", [](OpSum const &ops) { return to_string(ops); });
}
} // namespace xdiag::julia
