#include "symmetrize.hpp"

#include <xdiag/operators/logic/permute.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/utils/scalar.hpp>

namespace xdiag {

OpSum symmetrize(Op const &op, PermutationGroup const &group) try {
  auto ops = OpSum({op});
  return symmetrize(ops, group);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum symmetrize(Op const &op, Representation const &irrep) try {
  auto ops = OpSum({op});
  return symmetrize(ops, irrep);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum symmetrize(OpSum const &ops, PermutationGroup const &group) try {
  return symmetrize(ops, Representation(group));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename T>
static OpSum symmetrize(OpSum const &ops, PermutationGroup const &group,
                        arma::Col<T> const &characters) try {
  OpSum ops_sym;
  int64_t N_group = group.size();
  for (auto [cpl, op] : ops.plain()) {
    // Create all symmetrized ops
    for (int64_t i = 0; i < N_group; ++i) {
      Scalar cpls = cpl.scalar();
      T bloch = characters(i);
      Op op_perm = permute(op, group[i]);
      Scalar cpl_sym = cpls * bloch / (T)N_group;
      ops_sym += cpl_sym * op_perm;
    }
  }
  return ops_sym;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum symmetrize(OpSum const &ops, Representation const &irrep) try {
  if (isreal(irrep) && isreal(ops)) {
    auto characters = irrep.characters().as<arma::vec>();
    return symmetrize(ops, irrep.group(), characters);
  } else {
    auto characters = irrep.characters().as<arma::cx_vec>();
    return symmetrize(ops, irrep.group(), characters);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
