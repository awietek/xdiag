#include "imaginary_time_evolve.hpp"

#include <xdiag/algorithms/time_evolution/evolve_lanczos.hpp>
namespace xdiag {

State imaginary_time_evolve(OpSum const &H, State psi, double time,
                            double precision, double shift) try {
  imaginary_time_evolve_inplace(H, psi, time, precision, shift);
  return psi;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void imaginary_time_evolve_inplace(OpSum const &H, State &psi, double time,
                                   double precision, double shift) try {
  // minus sign in exp(-Ht) implemented here
  evolve_lanczos_inplace(H, psi, -time, precision, shift);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
