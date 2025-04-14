// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_USE_MPI
#include <xdiag/basis/spinhalf_distributed/basis_sz.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <class bit_t, typename coeff_t>
void transpose(BasisSz<bit_t> const &basis, coeff_t const *vec_in,
               bool reverse = true);

} // namespace xdiag::basis::spinhalf_distributed
#endif
