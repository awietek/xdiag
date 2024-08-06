#include "state.hpp"

namespace xdiag::julia {
void define_gpwf(jlcxx::Module &mod) {

  mod.add_type<GPWF>("cxx_GPWF")
      .constructor<>()
      .constructor<arma::mat const &, int64_t>()
      .constructor<arma::cx_mat const &, int64_t>()
      .method("n_sites",
              [](GPWF const &s) { JULIA_XDIAG_CALL_RETURN(s.n_sites()) })
      .method("n_up", [](GPWF const &s) { JULIA_XDIAG_CALL_RETURN(s.n_up()) })
      .method("isreal",
              [](GPWF const &s) { JULIA_XDIAG_CALL_RETURN(s.isreal()) });

  mod.method("to_string", [](GPWF const &s) { return to_string(s); });
}

} // namespace xdiag::julia
