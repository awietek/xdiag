#pragma once

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::operators {

OpSum clean_zeros(OpSum const &ops, double precision = 1e-12);

// Checks for the correct format of ops
void check_op(Op const &op, int64_t n_sites_total, int64_t n_sites_op,
              bool disjoint, std::string type);

void check_op_in_range(Op const &op, int64_t n_sites);
void check_op_has_correct_number_of_sites(Op const &op, int64_t ns);
void check_op_has_disjoint_sites(Op const &op);
void check_op_coupling_has_type(Op const &op, std::string type);
void check_op_coupling_has_type(Op const &op, std::string type1,
                                std::string type2);

} // namespace xdiag::operators
