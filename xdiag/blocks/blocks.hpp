// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <cstdint>
#include <ostream>
#include <string>
#include <variant>

#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {
using Block = std::variant<Spinhalf, Boson, Fermion>;

XDIAG_API int64_t dim(Block const &block);
XDIAG_API int64_t size(Block const &block);
XDIAG_API int64_t nsites(Block const &block);
XDIAG_API bool isreal(Block const &block);

// Output block obtained by applying ops to block_in: the quantum-number sectors
// of block_in are shifted by the representation ops transforms under (U(1)
// charges add, permutation characters multiply). Throws if ops has no
// well-defined sector under one of the block's symmetries.
Spinhalf block(OpSum const &ops, Spinhalf const &block_in);
Boson block(OpSum const &ops, Boson const &block_in);
Fermion block(OpSum const &ops, Fermion const &block_in);
Block block(OpSum const &ops, Block const &block_in);

// True if applying ops to block_in lands exactly in block_out, i.e. block_out
// is the output block inferred from the operator quantum numbers.
bool blocks_match(OpSum const &ops, Block const &block_in,
                  Block const &block_out);

XDIAG_API std::ostream &operator<<(std::ostream &out, Block const &block);
XDIAG_API std::string to_string(Block const &block);

// Render a ProductState using the labels appropriate to its block (spin-1/2
// arrows, boson occupation numbers on a blue->red gradient).
XDIAG_API std::string to_string(ProductState const &state, Block const &block);

} // namespace xdiag
