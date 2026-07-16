// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

// Validates that every Op type appearing in ops is contained in the algebra's
// allowed_types set. Throws if an unsupported type is encountered. An empty
// allowed_types set disables the check.
void check_allowed_types(OpSum const &ops, Algebra const &algebra);

// Validates that every site appearing in ops lies in [0, algebra.nsites).
// Site-less ops (e.g. Id, HubbardU) are skipped. Throws on an out-of-range
// site.
void check_sites_in_range(OpSum const &ops, Algebra const &algebra);

} // namespace xdiag::algebra
