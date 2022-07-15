#include "operator_qns.h"

#include <hydra/utils/logger.h>

namespace hydra::utils {

int spinhalf_nup(BondList const &bonds, Couplings const &cpls) {

  auto clean_bonds = clean_bondlist(bonds, cpls, spinhalf_bonds);
  if (clean_bonds.size() == 0)
    return 0;

  auto type = clean_bonds[0].type();
  if (spinhalf_bond_nup.find(type) == spinhalf_bond_nup.end())
    Log.warn("Warning: unknown bond type ignored: {}!", type);

  int nup0 = spinhalf_bond_nup.at(clean_bonds[0].type());
  for (auto bond : clean_bonds) {
    auto type = bond.type();
    if (spinhalf_bond_nup.find(type) == spinhalf_bond_nup.end())
      Log.warn("Warning: unknown bond type ignored: {}!", type);

    if (spinhalf_bond_nup.at(type) != nup0)
      return undefined_qn;
  }
  return nup0;
}

std::pair<int, int> electron_nup_ndn(BondList const &bonds,
                                              Couplings const &cpls) {

  auto clean_bonds = clean_bondlist(bonds, cpls, electron_bonds);
  if (clean_bonds.size() == 0)
    return {0, 0};

  auto type = clean_bonds[0].type();
  if (electron_bond_nup_ndn.find(type) == electron_bond_nup_ndn.end())
    Log.err("Warning: unknown bond type ignored: {}!", type);

  auto nup_ndn0 = electron_bond_nup_ndn.at(clean_bonds[0].type());
  for (auto bond : clean_bonds) {
    auto type = bond.type();
    if (electron_bond_nup_ndn.find(type) == electron_bond_nup_ndn.end())
      Log.err("Warning: unknown bond type ignored: {}!", type);

    if (electron_bond_nup_ndn.at(type) != nup_ndn0)
      return {undefined_qn, undefined_qn};
  }
  return nup_ndn0;
}

std::pair<int, int> tj_nup_ndn(BondList const &bonds,
                                        Couplings const &cpls) {
  return electron_nup_ndn(bonds, cpls);
}

} // namespace hydra::utils
