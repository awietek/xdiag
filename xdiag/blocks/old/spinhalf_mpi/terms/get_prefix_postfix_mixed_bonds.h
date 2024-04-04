#pragma once

#include <algorithm>

#include <xdiag/operators/bondlist.h>
#include <xdiag/operators/couplings.h>

namespace xdiag::terms {

std::tuple<BondList, BondList, BondList>
get_prefix_postfix_mixed_bonds(BondList const &bonds, int n_postfix_sites);

} // namespace xdiag::terms
