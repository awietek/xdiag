#pragma once

#include <utility>

#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra::terms::spinhalf {

const std::vector<std::string> special_bond_types = {
    "HB", "HEISENBERG", "EXCHANGE", "ISING",
    "SZ", "S+",         "S-",       "SCALARCHIRALITY"};

BondList compile(BondList const &bonds, double precision = 1e-12);

} // namespace hydra::terms::spinhalf
