#pragma once
#include <string>
#include <vector>

#include <xdiag/operators/bondlist.hpp>

namespace xdiag::tj {

BondList compile(BondList const &bonds, int64_t n_sites, double precision);

} // namespace xdiag::tj
