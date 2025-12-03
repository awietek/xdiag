// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace xdiag::symmetries {

template <typename bit_t, typename int_t>
bit_t apply_permutation(bit_t state, int_t nsites, const int_t *permutation);

}
