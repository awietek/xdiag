#include "electron_symmetric_v2.h"

namespace hydra {

template <class bit_t, class GroupAction>
std::pair<std::vector<int>, std::vector<std::pair<idx_t, idx_t>>>
compute_stabilizer_arrays(IndexingSymmetricFermionic const &idxes,
                          GroupAction const &group_action) {

  std::vector<int> stabilizers;
  std::vector<std::pair<idx_t, idx_t>> stabilizers_limits;

  for (idx_t idx = 0; idx < idxes.size(); ++idx) {
    bit_t rep = idxes.state(idx);

    std::vector<int> stable_syms;
    for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
      bit_t trans = group_action.apply(sym, rep);
      if (trans == rep)
        stable_syms.push_back(sym);
    }

    idx_t begin = idxes_up_.size();
    idx_t end = begin + stable_syms.size();
    stabilizers.insert(stabilizers_up_.end(), stable_syms.begin(),
                       stable_syms.end());
    stabilizers_limits.push_back({begin, end});
  }

  return {stabilizers, stabilizers_limits};
}

template <class bit_t, class GroupAction>
ElectronSymmetricV2<bit_t, GroupAction>::ElectronSymmetricV2(
    int n_sites, int nup, int ndn, PermutationGroup permutation_group,
    Representation irrep)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      permutation_group_(
          (irrep.allowed_symmetries().size() > 0)
              ? permutation_group.subgroup(irrep.allowed_symmetries())
              : permutation_group),
      group_action_(permutation_group_), irrep_(irrep),
      idxes_up_(n_sites, nup, permutation_group_, irrep),
      idxes_dn_(n_sites, ndn, permutation_group_, irrep) {
  
  std::tie(stabilizers_up_, stabilizers_up_limits_) =
      compute_stabilizer_arrays(idxes_up_, group_action_);
  std::tie(stabilizers_dn_, stabilizers_dn_limits_) =
      compute_stabilizer_arrays(idxes_dn_, group_action_);
}

} // namespace hydra
