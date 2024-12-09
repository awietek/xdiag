#include "symmetrize.hpp"

namespace xdiag {

OpSum symmetrize(Op const &op, PermutationGroup const &group) {
  auto ops = OpSum({op});
  return symmetrize(ops, group);
}

OpSum symmetrize(Op const &op, PermutationGroup const &group,
                 Representation const &irrep) {
  auto ops = OpSum({op});
  return symmetrize(ops, group, irrep);
}

OpSum symmetrize(OpSum const &ops, PermutationGroup const &group) {
  return symmetrize(ops, group, trivial_representation(group.size()));
}

// TODO: needs to be adapted for MATRIX type
OpSum symmetrize(OpSum const &ops, PermutationGroup const &group,
                 Representation const &irrep) try {
  OpSum ops_sym;
  int64_t N_group = group.size();
  for (auto [cpl, op] : ops.plain()) {

    std::string type = op.type();

    // Create all symmetrized ops
    for (int64_t i = 0; i < N_group; ++i) {
      Permutation perm = group[i];
      complex bloch = irrep.character(i);

      std::vector<int64_t> sites_sym(op.size(), 0);
      for (int64_t site_idx = 0; site_idx < op.size(); ++site_idx) {
        sites_sym[site_idx] = perm[op[site_idx]];
      }
      Scalar cpl_sym(bloch * cpl.scalar().as<complex>() / (complex)N_group);
      ops_sym += cpl_sym * Op(type, sites_sym);
    }
  }
  return ops_sym;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
