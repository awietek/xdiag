#include "op.hpp"

namespace xdiag::julia {

void define_op(jlcxx::Module &mod) {

  mod.add_type<Op>("cxx_Op")
      .constructor<>()
      .constructor<std::string>()

      .constructor<std::string, int64_t>()
      .constructor<std::string, std::vector<int64_t> const &>()

      .constructor<std::string, int64_t, arma::mat const &>()
      .constructor<std::string, std::vector<int64_t> const &,
                   arma::mat const &>()

      .constructor<std::string, int64_t, arma::cx_mat const &>()
      .constructor<std::string, std::vector<int64_t> const &,
                   arma::cx_mat const &>()

      .method("type", [](Op const &op) { JULIA_XDIAG_CALL_RETURN(op.type()) })
      .method("size", [](Op const &op) { JULIA_XDIAG_CALL_RETURN(op.size()) })
      .method("getindex",
              [](Op const &op, int64_t idx) {
                JULIA_XDIAG_CALL_RETURN(op[idx - 1])
              })
      .method("sites", [](Op const &op) { JULIA_XDIAG_CALL_RETURN(op.sites()) })
      .method("isreal",
              [](Op const &op) { JULIA_XDIAG_CALL_RETURN(isreal(op)) });

  mod.method("to_string", [](Op const &op) { return to_string(op); });
}

} // namespace xdiag::julia
