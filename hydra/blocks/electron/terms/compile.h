#pragma once

#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra::electron {

const std::vector<std::string> special_bond_types = {
    "HB", "HEISENBERG", "EXCHANGE", "ISING", "HOP", "HOPUP", "HOPDN"};

BondList compile(BondList const &bonds, double precision = 1e-12);

} // namespace hydra::electron
