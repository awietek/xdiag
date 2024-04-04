#pragma once

#include <xdiag/common.h>
#include <xdiag/operators/bondlist.h>

namespace xdiag::spinhalf {

const std::vector<std::string> special_bond_types = {
    "HB", "HEISENBERG", "EXCHANGE", "ISING",
    "SZ", "S+",         "S-",       "SCALARCHIRALITY"};

BondList compile(BondList const &bonds, double precision = 1e-12);

} // namespace xdiag::spinhalf
