#pragma once

#include <xdiag/common.hpp>
#include <xdiag/operators/bondlist.hpp>

namespace xdiag::electron {

const std::vector<std::string> special_bond_types = {
    "HB",       "EXCHANGE", "ISING", "HOP",    "HOPUP",  "HOPDN", "NUMBER",
    "NUMBERUP", "NUMBERDN", "SZ",    "CDAGUP", "CDAGDN", "CUP",   "CDN"};

BondList compile(BondList const &bonds, double precision = 1e-12);

} // namespace xdiag::electron
