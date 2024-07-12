#pragma once
#include <string>
#include <vector>

#include <xdiag/operators/bondlist.hpp>

namespace xdiag::tj {

BondList compile(BondList const &bonds, double precision = 1e-12);

} // namespace xdiag::tj
