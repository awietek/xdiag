#include "qns.hpp"

#include <algorithm>

#include <xdiag/operators/compiler.hpp>
#include <xdiag/operators/non_branching_bonds.hpp>

namespace xdiag::spinhalf {

int64_t nup(BondList bonds, double precision) {
  bonds = operators::compile_explicit(bonds, precision, "keep");
  bonds = operators::non_branching_bonds(bonds, precision);

  int64_t nup = 0;
  bool first_bond = true;
  for (Bond bond : bonds) {

    int64_t nup_bond = 0;
    if (bond.type_defined()) {
      std::string type = bond.type();
      if (special_bonds_nup.count(type)) {
        nup_bond = special_bonds_nup.at(type);
      } else {
        Log.err("Error: cannot determine nup of bondlist.");
      }
    } else {
      auto bond_nb =
          operators::NonBranchingBond<uint64_t, complex>(bond, precision);
      nup_bond = bond_nb.number_difference();
    }
    if (first_bond) {
      nup = nup_bond;
      first_bond = false;
    } else {
      if (nup_bond != nup) {
        Log("xxx");
        return undefined_qn;
      }
    }
  }
  return nup;
}
} // namespace xdiag::spinhalf
