// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "blocks.hpp"

#include <type_traits>

#include <xdiag/algebra/algebras/spinhalf_implementation_algebra.hpp>
#include <xdiag/algebra/representation.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>

namespace xdiag {

// Symmetry sectors of the output block obtained by applying ops to a Spinhalf:
// the input sectors shifted by the representation ops transforms under (U(1)
// charges add, permutation characters multiply). Does not build a block.
static RepresentationSet output_irreps(OpSum const &ops,
                                       Spinhalf const &block_in) try {
  // How ops transforms under each symmetry sector of the input block. A sector
  // is omitted from `shift` whenever ops has no well-defined quantum number for
  // it.
  RepresentationSet in = block_in.irreps();
  RepresentationSet shift = algebra::representations(
      ops, in, algebra::spinhalf_implementation_algebra());

  // Every symmetry of the input block must be matched, otherwise applying ops
  // does not map block_in into a single output block.
  for (Representation const &rep : in) {
    if (!shift.has_type(rep.type())) {
      XDIAG_THROW(fmt::format(
          "Cannot determine output block: the OpSum has no well-defined "
          "quantum number for the symmetry of type \"{}\" carried by the input "
          "block.",
          rep.type()));
    }
  }
  return in * shift;
}
XDIAG_CATCH

Spinhalf block(OpSum const &ops, Spinhalf const &block_in) try {
  RepresentationSet out = output_irreps(ops, block_in);

  // If ops leaves every sector unchanged the output block equals the input
  // block; return a shallow copy (shares the basis) instead of rebuilding it.
  // isapprox, not ==, because the computed characters only match approximately.
  if (block_in.irreps().isapprox(out)) {
    return block_in;
  }
  return Spinhalf(block_in.nsites(), out);
}
XDIAG_CATCH

Block block(OpSum const &ops, Block const &block_in) try {
  return std::visit([&](auto const &b) -> Block { return block(ops, b); },
                    block_in);
}
XDIAG_CATCH

bool blocks_match(OpSum const &ops, Block const &block_in,
                  Block const &block_out) try {
  // Check that block_out is the sector ops maps block_in into, comparing only
  // nsites and irreps -- no block is copied or constructed. isapprox, not ==,
  // because the output characters are computed and only match approximately.
  return std::visit(
      [&](auto const &bin, auto const &bout) -> bool {
        if constexpr (std::is_same_v<decltype(bin), decltype(bout)>) {
          return (bin.nsites() == bout.nsites()) &&
                 isapprox(output_irreps(ops, bin), bout.irreps());
        } else { // different block types never match
          return false;
        }
      },
      block_in, block_out);
}
XDIAG_CATCH

int64_t dim(Block const &block) {
  return std::visit([&](auto &&b) { return b.dim(); }, block);
}

int64_t size(Block const &block) {
  return std::visit([&](auto &&b) { return b.size(); }, block);
}

int64_t nsites(Block const &block) {
  return std::visit([&](auto &&b) { return b.nsites(); }, block);
}

bool isreal(Block const &block) {
  return std::visit([&](auto &&b) { return b.isreal(); }, block);
}

std::ostream &operator<<(std::ostream &out, Block const &block) {
  std::visit([&](auto &&block) { out << block; }, block);
  return out;
}
std::string to_string(Block const &block) { return to_string_generic(block); }

} // namespace xdiag
