#pragma once

#include <mpi.hpp>
#include <xdiag/operators/bondlist.hpp>
#include <xdiag/operators/couplings.hpp>

namespace xdiag {

std::tuple<int, int, double>
get_exchange_s1_s2_Jhalf_mpi(Bond const &bond, Couplings const &couplings);

}
