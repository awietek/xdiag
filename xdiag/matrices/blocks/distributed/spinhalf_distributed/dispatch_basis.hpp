// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_DISTRIBUTED

#include <cstdint>
#include <utility>

#include <xdiag/basis/distributed/basis_spinhalf_distributed.hpp>
#include <xdiag/blocks/distributed/spinhalf_distributed.hpp>
#include <xdiag/matrices/dispatcher.hpp>

namespace xdiag::matrices {

// Resolve the type-erased basis of a SpinhalfDistributed block to the concrete
// BasisSpinhalfDistributed<bit_t> backend and invoke f(basis_in, basis_out).
template <typename func_t>
void dispatch_basis(SpinhalfDistributed const &block_in,
                    SpinhalfDistributed const &block_out, func_t &&f) {
  dispatch_basis_types<basis::BasisSpinhalfDistributed<uint32_t>,
                       basis::BasisSpinhalfDistributed<uint64_t>>(
      *block_in.basis(), *block_out.basis(), std::forward<func_t>(f));
}

} // namespace xdiag::matrices
#endif
