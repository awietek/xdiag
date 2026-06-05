// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "convert_to_matrix_rule.hpp"

#include <vector>

#include <xdiag/algebra/utils/op_to_matrix_op.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

MonomialRule
convert_to_matrix_rule(std::set<std::string> const &protected_single_op_types) {
  return [protected_single_op_types](
             Monomial const &mono) -> std::optional<OpSum> {
    for (int64_t k = 0; k < mono.size(); ++k) {
      std::string const &type = mono[k].type();
      if (type == "Matrix") {
        continue;
      }
      if (type == "Id") {
        continue; // handled by id_absorption_rule
      }
      // In a size-1 monomial, protected types stay as-is.
      if (mono.size() == 1 && protected_single_op_types.count(type) > 0) {
        return std::nullopt;
      }
      Op mat_op = op_to_matrix_op(mono[k]);
      std::vector<Op> ops = mono.ops();
      ops[k] = mat_op;
      return OpSum(Monomial(ops));
    }
    return std::nullopt;
  };
}

} // namespace xdiag::algebra
