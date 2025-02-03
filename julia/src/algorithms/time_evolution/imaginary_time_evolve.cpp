#include "imaginary_time_evolve.hpp"

namespace xdiag::julia {

void define_imaginary_time_evolve(jlcxx::Module &mod) {

  mod.method("cxx_imaginary_time_evolve",
             [](OpSum const &ops, State psi, double time, double precision,
                double shift) {
               JULIA_XDIAG_CALL_RETURN(
                   imaginary_time_evolve(ops, psi, time, precision, shift));
             });

  mod.method("cxx_imaginarytime_evolve_inplace",
             [](OpSum const &ops, State &psi, double time, double precision,
                double shift) {
               JULIA_XDIAG_CALL_RETURN(imaginary_time_evolve_inplace(
                   ops, psi, time, precision, shift));
             });
}
} // namespace xdiag::julia
