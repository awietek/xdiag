#pragma once

#include <utility>

#include <xdiag/operators/bond.hpp>
#include <xdiag/operators/bondlist.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

BondList symmetrized_operator(Bond const &bond, PermutationGroup const &group);

BondList symmetrized_operator(Bond const &bond, PermutationGroup const &group,
                              Representation const &irrep);

BondList symmetrized_operator(BondList const &bonds,
                              PermutationGroup const &group);

BondList symmetrized_operator(BondList const &bonds,
                              PermutationGroup const &group,
                              Representation const &irrep);

} // namespace xdiag
