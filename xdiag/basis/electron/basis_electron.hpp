// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <variant>

#include <xdiag/basis/electron/basis_no_np.hpp>
#include <xdiag/basis/electron/basis_np.hpp>
#include <xdiag/basis/electron/basis_symmetric_no_np.hpp>
#include <xdiag/basis/electron/basis_symmetric_np.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis {

// clang-format off
using BasisElectron =
  std::variant<electron::BasisNp<uint32_t>,
	       electron::BasisNoNp<uint32_t>,
	       electron::BasisSymmetricNp<uint32_t>,
	       electron::BasisSymmetricNoNp<uint32_t>,
	       electron::BasisNp<uint64_t>,
	       electron::BasisNoNp<uint64_t>,
	       electron::BasisSymmetricNp<uint64_t>,
	       electron::BasisSymmetricNoNp<uint64_t>>;
// clang-format on

// clang-format off
using BasisElectronIterator =
  std::variant<electron::BasisNpIterator<uint32_t>,
	       electron::BasisNoNpIterator<uint32_t>,
	       electron::BasisSymmetricNpIterator<uint32_t>,
	       electron::BasisSymmetricNoNpIterator<uint32_t>,
	       electron::BasisNpIterator<uint64_t>,
	       electron::BasisNoNpIterator<uint64_t>,
	       electron::BasisSymmetricNpIterator<uint64_t>,
	       electron::BasisSymmetricNoNpIterator<uint64_t>>;
// clang-format on

int64_t dim(BasisElectron const &basis);
int64_t size(BasisElectron const &basis);
template <typename bit_t> bool has_bit_t(BasisElectron const &);

} // namespace xdiag::basis
