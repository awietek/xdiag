#include "get_prefix_postfix_mixed_bonds.h"

namespace hydra::terms {

std::tuple<BondList, BondList, BondList>
get_prefix_postfix_mixed_bonds(BondList const &bonds, int n_postfix_sites) {

  BondList prefix_bonds;
  BondList postfix_bonds;
  BondList mixed_bonds;
  for (auto bond : bonds) {
    auto sites = bond.sites();
    // Postfix bond
    if (std::all_of(sites.begin(), sites.end(), [&n_postfix_sites](int i) {
          return i < n_postfix_sites;
        })) {
      postfix_bonds << bond;
    }
    // Prefix bond
    else if (std::all_of(sites.begin(), sites.end(), [&n_postfix_sites](int i) {
               return i >= n_postfix_sites;
             })) {
      prefix_bonds << bond;
    }
    // Mixed bond
    else {
      mixed_bonds << bond;
    }
  }
  return {prefix_bonds, postfix_bonds, mixed_bonds};
}

} // namespace hydra::terms
