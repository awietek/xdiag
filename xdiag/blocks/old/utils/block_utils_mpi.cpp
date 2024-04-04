#include "block_utils_mpi.hpp"

#include <xdiag/mpi/logger_mpi.hpp>

namespace xdiag {

std::tuple<int, int, double>
get_exchange_s1_s2_Jhalf_mpi(Bond const &bond, Couplings const &couplings) {
  if (bond.size() != 2)
    LogMPI.err("Error computing Exchange: "
               "bond must have exactly two sites defined");
  double Jhalf = lila::real(couplings[bond.coupling()]) / 2.;
  int s1 = bond.site(0);
  int s2 = bond.site(1);
  if (s1 == s2)
    LogMPI.err("Error computing Exchange: "
               "operator acting on twice the same site");
  return {std::min(s1, s2), std::max(s1, s2), Jhalf};
}

} // namespace xdiag
