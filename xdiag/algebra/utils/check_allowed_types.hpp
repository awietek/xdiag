// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/algebras/algebra.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

// Validates that every Op type appearing in ops is contained in the algebra's
// allowed_types set. Throws if an unsupported type is encountered. An empty
// allowed_types set disables the check.
void check_allowed_types(OpSum const &ops, Algebra const &algebra);

} // namespace xdiag::algebra
