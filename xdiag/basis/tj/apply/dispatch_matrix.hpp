// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/blocks/tj.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis {

template <typename coeff_t>
void dispatch_matrix(OpSum const &ops, tJ const &block_in, tJ const &block_out,
                     coeff_t *mat, int64_t m);

} // namespace xdiag::basis
