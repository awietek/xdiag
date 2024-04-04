#pragma once

#include <mpi.h>
#include <xdiag/operators/bondlist.h>
#include <xdiag/operators/couplings.h>

namespace xdiag {

std::tuple<int, int, double>
get_exchange_s1_s2_Jhalf_mpi(Bond const &bond, Couplings const &couplings);

}
