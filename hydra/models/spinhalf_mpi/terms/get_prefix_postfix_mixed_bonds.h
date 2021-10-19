#pragma once

#include <algorithm>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms::spinhalf_mpi {

std::tuple<BondList, BondList, BondList>
get_prefix_postfix_mixed_bonds(BondList const &bonds, int n_postfix_sites);

} // namespace hydra::terms::spinhalf_mpi
