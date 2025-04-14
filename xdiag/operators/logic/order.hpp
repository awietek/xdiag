// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/scalar.hpp>

namespace xdiag {

std::pair<Scalar, Op> order(Op const &op);
std::pair<Scalar, Op> order(Scalar const &alpha, Op const &op);
OpSum order(OpSum const &ops);

bool less(Op const &o1, Op const &o2);
bool less(Scalar const &s1, Op const &s2);
bool less(Matrix const &m1, Matrix const &m2);

} // namespace xdiag
