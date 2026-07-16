// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_DISTRIBUTED

#include <cstdint>
#include <utility>

#include <xdiag/basis/distributed/basis_electron_distributed.hpp>
#include <xdiag/blocks/distributed/electron_distributed.hpp>
#include <xdiag/kernels/dispatcher.hpp>

namespace xdiag::kernels {

// Resolve the type-erased basis of a ElectronDistributed block to the concrete
// BasisElectronDistributed<bit_t> backend and invoke f(basis_in, basis_out).
template <typename func_t>
void dispatch_basis(ElectronDistributed const &block_in, ElectronDistributed const &block_out,
                    func_t &&f) {
  dispatch_basis_types<basis::BasisElectronDistributed<uint32_t>,
                       basis::BasisElectronDistributed<uint64_t>>(
      *block_in.basis(), *block_out.basis(), std::forward<func_t>(f));
}

} // namespace xdiag::kernels
#endif
