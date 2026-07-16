// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "matrix_rules.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/utils/combine_matrix_ops.hpp>
#include <xdiag/algebra/utils/op_to_matrix_op.hpp>
#include <xdiag/algebra/utils/permute_matrix_op.hpp>
#include <xdiag/algebra/utils/replace_pair.hpp>

namespace xdiag::algebra {

std::optional<OpSum> boson_expand(Op const &op, Algebra const &alg) {
  std::string t = op.type();
  if (t == "Hop") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Adag", i) * Op("A", j));
    r += -1.0 * (Op("Adag", j) * Op("A", i));
    return r;
  }
  if (t == "HopAsym") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Adag", i) * Op("A", j));
    r += +1.0 * (Op("Adag", j) * Op("A", i));
    return r;
  }
  if (t == "HubbardU") {
    // (1/2) sum_i N{i} (N{i} - 1)
    OpSum r;
    for (int64_t i = 0; i < alg.nsites; ++i) {
      r += 0.5 * (Op("N", i) * Op("N", i));
      r += -0.5 * Op("N", i);
    }
    return r;
  }
  if (t == "TotalN") {
    OpSum r;
    for (int64_t i = 0; i < alg.nsites; ++i) {
      r += Op("N", i);
    }
    return r;
  }
  return std::nullopt;
}

std::optional<OpSum> spinhalf_expand(Op const &op, Algebra const &alg) {
  std::string t = op.type();
  if (t == "SdotS") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Exchange", {i, j})) + OpSum(Op("SzSz", {i, j}));
  }
  if (t == "TotalSz") {
    OpSum r;
    for (int64_t i = 0; i < alg.nsites; ++i) {
      r += Op("Sz", i);
    }
    return r;
  }
  return std::nullopt;
}

std::optional<OpSum> matrix_simplify(Monomial const &mono, Algebra const &alg) {
  int64_t n = mono.size();

  // 1. Convert the first non-Matrix operator to a local Matrix operator. A sole
  //    operator whose type is kept named (it has a dedicated kernel) is left.
  for (int64_t k = 0; k < n; ++k) {
    std::string const &type = mono[k].type();
    if (type == "Matrix" || type == "Id") {
      continue;
    }
    if (n == 1 && alg.kept_named.count(type) > 0) {
      return std::nullopt;
    }
    std::vector<Op> ops = mono.ops();
    ops[k] = op_to_matrix_op(mono[k], alg.d);
    return OpSum(Monomial(ops));
  }

  // 2. Fuse two adjacent Matrix operators.
  for (int64_t k = 0; k + 1 < n; ++k) {
    if (mono[k].type() == "Matrix" && mono[k].hassites() &&
        mono[k + 1].type() == "Matrix" && mono[k + 1].hassites()) {
      Op combined = combine_matrix_ops({mono[k], mono[k + 1]}, alg.d);
      return replace_pair(mono, k, OpSum(combined));
    }
  }

  // 3. Sort the sites within a Matrix operator.
  for (int64_t k = 0; k < n; ++k) {
    Op const &op = mono[k];
    if (op.type() != "Matrix" || !op.hassites()) {
      continue;
    }
    std::vector<int64_t> const &sites = op.sites();
    if (std::is_sorted(sites.begin(), sites.end())) {
      continue;
    }
    std::vector<Op> ops = mono.ops();
    ops[k] = permute_matrix_op(op, alg.d);
    return OpSum(Monomial(ops));
  }
  return std::nullopt;
}

} // namespace xdiag::algebra
