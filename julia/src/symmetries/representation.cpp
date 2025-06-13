// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "representation.hpp"

namespace xdiag::julia {
void define_representation(jlcxx::Module &mod) {

  mod.add_type<Representation>("cxx_Representation")
      .method("isreal", &Representation::isreal)
      .method("size", &Representation::size);

  mod.method("multiply",
             [](Representation const &r1, Representation const &r2) {
               JULIA_XDIAG_CALL_RETURN(multiply(r1, r2));
             });

  mod.method("to_string", [](Representation const &r) { return to_string(r); });

  mod.method("construct_Representation",
             []() { JULIA_XDIAG_CALL_RETURN(Representation()); });
  mod.method("construct_Representation", [](PermutationGroup const &group) {
    JULIA_XDIAG_CALL_RETURN(Representation(group));
  });
  mod.method("construct_Representation",
             [](PermutationGroup const &group, arma::vec vec) {
               JULIA_XDIAG_CALL_RETURN(Representation(group, vec));
             });
  mod.method("construct_Representation",
             [](PermutationGroup const &group, arma::cx_vec vec) {
               JULIA_XDIAG_CALL_RETURN(Representation(group, vec));
             });
}

} // namespace xdiag::julia
