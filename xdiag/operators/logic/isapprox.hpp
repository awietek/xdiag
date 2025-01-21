#pragma once

#include <optional>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/scalar.hpp>

namespace xdiag {
bool isapprox(Op const &op1, Op const &op2, double rtol = 1e-12,
              double atol = 1e-12);
bool isapprox(OpSum const &ops1, OpSum const &ops2, double rtol = 1e-12,
              double atol = 1e-12);
std::optional<Scalar> isapprox_multiple(OpSum const &ops1, OpSum const &ops2,
                                        double rtol = 1e-12,
                                        double atol = 1e-12);
} // namespace xdiag
