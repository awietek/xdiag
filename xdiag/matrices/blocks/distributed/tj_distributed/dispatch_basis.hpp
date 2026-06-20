// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_DISTRIBUTED

#include <cstdint>
#include <utility>

#include <xdiag/basis/distributed/basis_tj_distributed.hpp>
#include <xdiag/blocks/distributed/tj_distributed.hpp>
#include <xdiag/matrices/dispatcher.hpp>

namespace xdiag::matrices {

// Resolve the type-erased basis of a tJDistributed block to the concrete
// BasistJDistributed<bit_t> backend and invoke f(basis_in, basis_out).
template <typename func_t>
void dispatch_basis(tJDistributed const &block_in, tJDistributed const &block_out,
                    func_t &&f) {
  dispatch_basis_types<basis::BasistJDistributed<uint32_t>,
                       basis::BasistJDistributed<uint64_t>>(
      *block_in.basis(), *block_out.basis(), std::forward<func_t>(f));
}

} // namespace xdiag::matrices
#endif
