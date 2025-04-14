// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_USE_MPI
#include <variant>
#include <xdiag/basis/electron_distributed/basis_np.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis {
// clang-format off
using BasisElectronDistributed =
  std::variant<electron_distributed::BasisNp<uint32_t>,
	       electron_distributed::BasisNp<uint64_t>>;
// clang-format on

// clang-format off
using BasisElectronDistributedIterator =
  std::variant<electron_distributed::BasisNpIterator<uint32_t>,
	       electron_distributed::BasisNpIterator<uint64_t>>;
// clang-format on

int64_t dim(BasisElectronDistributed const &basis);
int64_t size(BasisElectronDistributed const &basis);
int64_t size_max(BasisElectronDistributed const &basis);
int64_t size_min(BasisElectronDistributed const &basis);

template <typename bit_t> bool has_bit_t(BasisElectronDistributed const &);

} // namespace xdiag::basis

#endif
