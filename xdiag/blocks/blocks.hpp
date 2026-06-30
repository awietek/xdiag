// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <cstdint>
#include <ostream>
#include <string>
#include <variant>

#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#ifdef XDIAG_DISTRIBUTED
#include <xdiag/blocks/distributed/electron_distributed.hpp>
#include <xdiag/blocks/distributed/spinhalf_distributed.hpp>
#include <xdiag/blocks/distributed/tj_distributed.hpp>
#endif
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/variants.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

using BlockShared = std::variant<Spinhalf, Boson, Fermion, Electron, tJ>;

#ifdef XDIAG_DISTRIBUTED
using BlockDistributed =
    std::variant<SpinhalfDistributed, tJDistributed, ElectronDistributed>;
using Block = utils::variant_union<BlockShared, BlockDistributed>::type;
#else
using Block = BlockShared;
#endif

XDIAG_API int64_t dim(Block const &block);
XDIAG_API int64_t size(Block const &block);
XDIAG_API int64_t nsites(Block const &block);
XDIAG_API bool isreal(Block const &block);

XDIAG_API std::ostream &operator<<(std::ostream &out, Block const &block);
XDIAG_API std::string to_string(Block const &block);

// Render a ProductState using the labels appropriate to its block (spin-1/2
// arrows, boson occupation numbers on a blue->red gradient).
XDIAG_API std::string to_string(ProductState const &state, Block const &block);

// Semantic equality: same block type, number of sites, local dimension, and
// (approximately) the same quantum-number sectors / irreps. Unlike operator==,
// which for some blocks compares basis-pointer identity, this compares the
// physical content, so two independently built blocks for the same sector
// compare equal.
XDIAG_API bool isapprox(Block const &b1, Block const &b2, double tol = 1e-12);

// Output block obtained by applying ops to block_in: the quantum-number sectors
// of block_in are shifted by the representation ops transforms under (U(1)
// charges add, permutation characters multiply). Throws if ops has no
// well-defined sector under one of the block's symmetries.
// Concrete-block overload, one body for every block type. Explicitly
// instantiated in blocks.cpp for each block, so the definition (and its heavy
// output_irreps / symmetry_algebra includes) stays out of this header.
template <typename block_t>
block_t block(OpSum const &ops, block_t const &block_in);
XDIAG_API Block block(OpSum const &ops, Block const &block_in);

// True if applying ops to block_in lands exactly in block_out, i.e. block_out
// is the output block inferred from the operator quantum numbers.
bool blocks_match(OpSum const &ops, Block const &block_in,
                  Block const &block_out);

template <typename T> struct is_distributed : std::false_type {};
#ifdef XDIAG_DISTRIBUTED
template <> struct is_distributed<SpinhalfDistributed> : std::true_type {};
template <> struct is_distributed<tJDistributed> : std::true_type {};
template <> struct is_distributed<ElectronDistributed> : std::true_type {};
#endif
template <typename T>
inline constexpr bool is_distributed_v = is_distributed<T>::value;

inline bool isdistributed(Block const &block) {
  return std::visit(
      [](auto &&blk) {
        using block_t = std::decay_t<decltype(blk)>;
        return is_distributed_v<block_t>;
      },
      block);
}

} // namespace xdiag
