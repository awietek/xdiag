// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/rewrite.hpp>

namespace xdiag::operators {

// MonomialRule: remove a site-free "Id" op from a monomial when other operators
// surround it. A standalone {Id} monomial (length 1) is left unchanged — it
// represents the scalar identity and must not be collapsed to an empty
// monomial.
MonomialRule id_absorption_rule();

} // namespace xdiag::operators
