#pragma once

#include <utility>

#include <xdiag/operators/bond.h>
#include <xdiag/operators/bondlist.h>
#include <xdiag/symmetries/permutation_group.h>
#include <xdiag/symmetries/representation.h>

namespace xdiag {

BondList symmetrized_operator(Bond const &bond,
                              PermutationGroup const &group);

BondList symmetrized_operator(Bond const &bond,
                              PermutationGroup const &group,
                              Representation const &irrep);
  
BondList symmetrized_operator(BondList const &bonds,
                              PermutationGroup const &group);

BondList symmetrized_operator(BondList const &bonds,
                              PermutationGroup const &group,
                              Representation const &irrep);

} // namespace xdiag
