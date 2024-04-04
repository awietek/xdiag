#pragma once

#include <xdiag/common.hpp>
#include <xdiag/operators/bondlist.hpp>

namespace xdiag::spinhalf {

const std::vector<std::string> special_bond_types = {
    "HB", "HEISENBERG", "EXCHANGE", "ISING",
    "SZ", "S+",         "S-",       "SCALARCHIRALITY"};

BondList compile(BondList const &bonds, double precision = 1e-12);

} // namespace xdiag::spinhalf
