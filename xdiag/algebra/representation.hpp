// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/symmetries/representation_set.hpp>

namespace xdiag::algebra {

std::optional<Representation> representation(OpSum const &ops,
                                             Representation const &irrep,
                                             Algebra const &algebra,
                                             double tol = 1e-12);

RepresentationSet representations(OpSum const &ops,
                                  RepresentationSet const &irreps,
                                  Algebra const &algebra, double tol = 1e-12);

} // namespace xdiag::algebra
