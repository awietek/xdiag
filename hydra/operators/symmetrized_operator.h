#pragma once

#include <utility>

#include <hydra/operators/bond.h>
#include <hydra/operators/bondlist.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra {

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

} // namespace hydra
