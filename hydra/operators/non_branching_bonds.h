#pragma once

#include <hydra/operators/bond.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

BondList non_branching_bonds(Bond const &bond, double precision = 1e-12);
BondList non_branching_bonds(BondList const &bonds, double precision = 1e-12);

} // namespace hydra
