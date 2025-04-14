// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "time_evolve.hpp"

namespace xdiag::julia {

void define_time_evolve(jlcxx::Module &mod) {

  mod.method("cxx_time_evolve", [](OpSum const &ops, State psi, double time,
                                   double precision, std::string algorithm) {
    JULIA_XDIAG_CALL_RETURN(time_evolve(ops, psi, time, precision, algorithm));
  });

  mod.method("cxx_time_evolve_inplace",
             [](OpSum const &ops, State &psi, double time, double precision,
                std::string algorithm) {
               JULIA_XDIAG_CALL_RETURN(
                   time_evolve_inplace(ops, psi, time, precision, algorithm));
             });
}
} // namespace xdiag::julia
