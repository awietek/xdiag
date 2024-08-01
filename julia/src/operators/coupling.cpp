#include "coupling.hpp"

namespace xdiag::julia {

void define_coupling(jlcxx::Module &mod) {
  mod.add_type<Coupling>("cxx_Coupling")
      .constructor<>()
      .constructor<std::string>()
      .constructor<double>()
      .constructor<complex>()
      .constructor<arma::mat>()
      .constructor<arma::cx_mat>()
      .method("isreal",
              [](Coupling const &cpl) { JULIA_XDIAG_CALL_RETURN(cpl.isreal()) })
      .method(
          "ismatrix",
          [](Coupling const &cpl) { JULIA_XDIAG_CALL_RETURN(cpl.ismatrix()) })
      .method("isexplicit",
              [](Coupling const &cpl) {
                JULIA_XDIAG_CALL_RETURN(cpl.isexplicit());
              })
      .method("type",
              [](Coupling const &cpl) { JULIA_XDIAG_CALL_RETURN(cpl.type()) })
      .method("is_string",
              [](Coupling const &cpl) {
                JULIA_XDIAG_CALL_RETURN(cpl.is<std::string>())
              })
      .method(
          "is_double",
          [](Coupling const &cpl) { JULIA_XDIAG_CALL_RETURN(cpl.is<double>()) })
      .method("is_complex",
              [](Coupling const &cpl) {
                JULIA_XDIAG_CALL_RETURN(cpl.is<complex>())
              })
      .method("is_mat",
              [](Coupling const &cpl) {
                JULIA_XDIAG_CALL_RETURN(cpl.is<arma::mat>())
              })
      .method("is_cx_mat",
              [](Coupling const &cpl) {
                JULIA_XDIAG_CALL_RETURN(cpl.is<arma::cx_mat>())
              })
      .method("as_string",
              [](Coupling const &cpl) {
                JULIA_XDIAG_CALL_RETURN(cpl.as<std::string>())
              })
      .method(
          "as_double",
          [](Coupling const &cpl) { JULIA_XDIAG_CALL_RETURN(cpl.as<double>()) })
      .method("as_complex",
              [](Coupling const &cpl) {
                JULIA_XDIAG_CALL_RETURN(cpl.as<complex>())
              })
      .method("as_mat",
              [](Coupling const &cpl) {
                JULIA_XDIAG_CALL_RETURN(cpl.as<arma::mat>())
              })
      .method("as_cx_mat", [](Coupling const &cpl) {
        JULIA_XDIAG_CALL_RETURN(cpl.as<arma::cx_mat>())
      });

  mod.method("to_string", [](Coupling const &cpl) {
    JULIA_XDIAG_CALL_RETURN(to_string(cpl))
  });
}

} // namespace xdiag::julia
