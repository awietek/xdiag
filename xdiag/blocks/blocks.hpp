// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <variant>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/common.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/blocks/electron_distributed.hpp>
#include <xdiag/blocks/spinhalf_distributed.hpp>
#include <xdiag/blocks/tj_distributed.hpp>
#endif

namespace xdiag {

#ifdef XDIAG_USE_MPI
using Block = std::variant<Spinhalf, tJ, Electron, SpinhalfDistributed,
                           tJDistributed, ElectronDistributed>;
#else
using Block = std::variant<Spinhalf, tJ, Electron>;
#endif

int64_t dim(Block const &block);
int64_t size(Block const &block);
int64_t nsites(Block const &block);
bool isreal(Block const &block);

template <typename block_t> constexpr bool isdistributed() {
  return !((std::is_same<block_t, Spinhalf>::value) ||
           (std::is_same<block_t, tJ>::value) ||
           (std::is_same<block_t, Electron>::value));
}

constexpr bool isdistributed(Block const &block) {
  return std::visit(
      [&](auto &&block) -> bool {
        using block_t = typename std::decay<decltype(block)>::type;
        return isdistributed<block_t>();
      },
      block);
}

std::ostream &operator<<(std::ostream &out, Block const &block);
std::string to_string(Block const &block);

} // namespace xdiag
