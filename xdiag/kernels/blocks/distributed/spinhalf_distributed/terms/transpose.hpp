// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_DISTRIBUTED
#include <xdiag/basis/distributed/basis_spinhalf_distributed.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <class bit_t, typename coeff_t>
void transpose(BasisSpinhalfDistributed<bit_t> const &basis, coeff_t const *vec_in,
               bool reverse = true);

} // namespace xdiag::basis::spinhalf_distributed
#endif
