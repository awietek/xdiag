#pragma once
#include <string>
#include <vector>

#include <hydra/operators/bondlist.h>

namespace hydra::tj {

const std::vector<std::string> special_bond_types = {
    "HB",     "HEISENBERG",   "EXCHANGE", "ISING",    "TJISING",
    "TJHB",   "TJHEISENBERG", "ISING",    "HOP",      "HOPUP",
    "HOPDN",  "NUMBER",       "NUMBERUP", "NUMBERDN", "SZ",
    "CDAGUP", "CDAGDN",       "CUP",      "CDN"};

BondList compile(BondList const &bonds, double precision = 1e-12);

} // namespace hydra::tj
