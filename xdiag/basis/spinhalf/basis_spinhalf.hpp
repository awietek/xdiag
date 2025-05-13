// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <variant>

#include <xdiag/basis/spinhalf/basis_no_sz.hpp>
#include <xdiag/basis/spinhalf/basis_sublattice.hpp>
#include <xdiag/basis/spinhalf/basis_symmetric_no_sz.hpp>
#include <xdiag/basis/spinhalf/basis_symmetric_sz.hpp>
#include <xdiag/basis/spinhalf/basis_sz.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis {

// clang-format off
using BasisSpinhalf = std::variant<
    spinhalf::BasisSz<uint32_t>,
    spinhalf::BasisNoSz<uint32_t>,
    spinhalf::BasisSymmetricSz<uint32_t>,
    spinhalf::BasisSymmetricNoSz<uint32_t>,
    spinhalf::BasisSublattice<uint32_t, 1>,
    spinhalf::BasisSublattice<uint32_t, 2>,
    spinhalf::BasisSublattice<uint32_t, 3>,
    spinhalf::BasisSublattice<uint32_t, 4>,
    spinhalf::BasisSublattice<uint32_t, 5>,
    spinhalf::BasisSz<uint64_t>,
    spinhalf::BasisNoSz<uint64_t>,
    spinhalf::BasisSymmetricSz<uint64_t>,
    spinhalf::BasisSymmetricNoSz<uint64_t>,
    spinhalf::BasisSublattice<uint64_t, 1>,
    spinhalf::BasisSublattice<uint64_t, 2>,
    spinhalf::BasisSublattice<uint64_t, 3>,
    spinhalf::BasisSublattice<uint64_t, 4>,
    spinhalf::BasisSublattice<uint64_t, 5>>;
// clang-format on

// clang-format off
using BasisSpinhalfIterator = std::variant<
  combinatorics::SubsetsIterator<uint32_t>,
  combinatorics::CombinationsIterator<uint32_t>,
  typename std::vector<uint32_t>::const_iterator,
  combinatorics::SubsetsIterator<uint64_t>,
  combinatorics::CombinationsIterator<uint64_t>,
  typename std::vector<uint64_t>::const_iterator>;
// clang-format on

int64_t dim(BasisSpinhalf const &basis);
int64_t size(BasisSpinhalf const &basis);
template <typename bit_t> bool has_bit_t(BasisSpinhalf const &);

} // namespace xdiag::basis
