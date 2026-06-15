// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <ostream>

#include <xdiag/extern/fmt/color.hpp>

namespace xdiag {

// Shared printer for all blocks (Spinhalf, Boson, Fermion, Electron). Given the
// output stream and the block, it renders the common, identically-formatted
// header and derives everything (name, header color, conserved charges, basis,
// hashes, dimension) from the block:
//
//   <block name>      (in the block's header color)
//   | nsites   : ...
//   | d        : ...
//   | <charge> : ...  (the block's conserved quantum numbers; "not conserved"
//   |                  when the charge is absent)
//   | permutation symmetries used   (only with a SitePermutation symmetry)
//   | irrep ID : ...                (only with a SitePermutation symmetry)
//   | basis    : ...  (basis type name with xdiag namespace prefixes stripped)
//   | ID       : ...  (block hash, hex)
//   | dimension: ...
//
// Defined in print_block.cpp and explicitly instantiated for the four blocks.
template <typename block_t>
void print_block(std::ostream &out, block_t const &block);

} // namespace xdiag
