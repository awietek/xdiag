#include "op.hpp"

namespace xdiag::julia {

void define_op(jlcxx::Module &mod) {
  mod.add_type<Op>("cxx_Op")
      .constructor<>()
      .constructor<std::string, Coupling const &,
                   std::vector<int64_t> const &>()
      .constructor<std::string, Coupling const &, int64_t>()
      .constructor<std::string, std::string, std::vector<int64_t> const &>()
      .constructor<std::string, std::string, int64_t>()
      .constructor<std::string, double, std::vector<int64_t> const &>()
      .constructor<std::string, double, int64_t>()
      .constructor<std::string, complex, std::vector<int64_t> const &>()
      .constructor<std::string, complex, int64_t>()
      .constructor<std::string, arma::mat const &,
                   std::vector<int64_t> const &>()
      .constructor<std::string, arma::mat const &, int64_t>()
      .constructor<std::string, arma::cx_mat const &,
                   std::vector<int64_t> const &>()
      .constructor<std::string, arma::cx_mat const &, int64_t>()
      .method("type", [](Op const &op) { JULIA_XDIAG_CALL_RETURN(op.type()) })
      .method("coupling",
              [](Op const &op) { JULIA_XDIAG_CALL_RETURN(op.coupling()) })
      .method("size", [](Op const &op) { JULIA_XDIAG_CALL_RETURN(op.size()) })
      .method("getindex", [](Op const &op,
                             int64_t idx) { JULIA_XDIAG_CALL_RETURN(op[idx-1]) })
      .method("sites", [](Op const &op) { JULIA_XDIAG_CALL_RETURN(op.sites()) })
      .method("isreal",
              [](Op const &op) { JULIA_XDIAG_CALL_RETURN(op.isreal()) })
      .method("ismatrix",
              [](Op const &op) { JULIA_XDIAG_CALL_RETURN(op.ismatrix()) })
      .method("isexplicit",
              [](Op const &op) { JULIA_XDIAG_CALL_RETURN(op.isexplicit()) });

  mod.method("to_string", [](Op const &op) { return to_string(op); });
}

} // namespace xdiag::julia
