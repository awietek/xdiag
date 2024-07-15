#pragma once

#include <xdiag/common.hpp>
#include <xdiag/operators/bondlist.hpp>

namespace xdiag::spinhalf {

BondList compile(BondList const &bonds, int64_t n_sites, double precision);

} // namespace xdiag::spinhalf
