#pragma once

#include <mpi.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

std::tuple<int, int, double>
get_exchange_s1_s2_Jhalf_mpi(Bond const &bond, Couplings const &couplings);

}
