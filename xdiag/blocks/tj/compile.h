#pragma once
#include <string>
#include <vector>

#include <xdiag/operators/bondlist.h>

namespace xdiag::tj {

const std::vector<std::string> special_bond_types = {
    "HB",     "HEISENBERG",   "EXCHANGE", "ISING",    "TJISING",
    "TJHB",   "TJHEISENBERG", "ISING",    "HOP",      "HOPUP",
    "HOPDN",  "NUMBER",       "NUMBERUP", "NUMBERDN", "SZ",
    "CDAGUP", "CDAGDN",       "CUP",      "CDN"};

BondList compile(BondList const &bonds, double precision = 1e-12);

} // namespace xdiag::tj
