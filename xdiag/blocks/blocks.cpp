// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "blocks.hpp"

#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
#include <type_traits>

#include <xdiag/algebra/algebras/algebra.hpp>
#include <xdiag/algebra/representation.hpp>
#include <xdiag/extern/fmt/color.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>
#include <xdiag/utils/variants.hpp>

namespace xdiag {

template <typename block_t>
static RepresentationSet output_irreps(OpSum const &ops,
                                       block_t const &block_in,
                                       algebra::Algebra const &algebra) try {
  RepresentationSet in = block_in.irreps();
  RepresentationSet shift = algebra::representations(ops, in, algebra);
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

template <typename block_t>
block_t block(OpSum const &ops, block_t const &block_in) try {
  RepresentationSet irreps_out =
      output_irreps(ops, block_in, algebra::symmetry_algebra(block_in));

  if (block_in.irreps().isapprox(irreps_out)) {
    return block_in;
  }
  if constexpr (std::is_same_v<block_t, Boson>) {
    // Boson carries the local dimension d explicitly in its constructor.
    return Boson(block_in.nsites(), block_in.d(), irreps_out);
  } else {
    return block_t(block_in.nsites(), irreps_out);
  }
}
XDIAG_CATCH

template Spinhalf block(OpSum const &, Spinhalf const &);
template Boson block(OpSum const &, Boson const &);
template Fermion block(OpSum const &, Fermion const &);
template Electron block(OpSum const &, Electron const &);

Block block(OpSum const &ops, Block const &block_in) try {
  // Generic visitor: the inner block(ops, b) resolves to the block<block_t>
  // template above. A block type whose template is not explicitly instantiated
  // (just above) is an undefined reference at link time.
  return std::visit([&](auto const &b) -> Block { return block(ops, b); },
                    block_in);
}
XDIAG_CATCH

bool blocks_match(OpSum const &ops, Block const &block_in,
                  Block const &block_out) try {
  return std::visit(
      [&](auto const &bin, auto const &bout) -> bool {
        if constexpr (std::is_same_v<decltype(bin), decltype(bout)>) {
          return isapprox(
              output_irreps(ops, bin, algebra::symmetry_algebra(block_in)),
              bout.irreps());
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

// Render a ProductState with the labels appropriate to its block. The block is
// what knows how to interpret the per-site integers: spin-1/2 maps them to
// colored arrows, bosons just print the occupation numbers. Highest site first,
// matching the bit ordering.
static std::string to_string(ProductState const &state, Spinhalf const &) {
  std::stringstream ss;
  for (int64_t i = state.size() - 1; i >= 0; --i) {
    if (state[i] == 1) { // Up
      ss << fmt::format(fg(fmt::color::light_blue), "↑");
    } else { // Dn
      ss << fmt::format(fg(fmt::color::orange), "↓");
    }
  }
  return ss.str();
}

// Bosons: occupation numbers in a fixed-width column (width set by the largest
// representable occupation d()-1, so e.g. a value of 5 lines up under 12), each
// colored on a blue->red gradient across the d() possible occupations.
static std::string to_string(ProductState const &state, Boson const &boson) {
  int64_t d = boson.d();
  int64_t max_occupation = (d > 0) ? d - 1 : 0;
  int width = static_cast<int>(std::to_string(max_occupation).size());
  std::stringstream ss;
  for (int64_t i = state.size() - 1; i >= 0; --i) {
    int64_t v = state[i];
    double f =
        (d > 1) ? static_cast<double>(v) / static_cast<double>(d - 1) : 0.0;
    uint8_t r = static_cast<uint8_t>(std::lround(255.0 * f));
    uint8_t b = static_cast<uint8_t>(std::lround(255.0 * (1.0 - f)));
    ss << fmt::format(fg(fmt::rgb(r, 0, b)), "{:>{}}", v, width);
    if (i > 0) {
      ss << " ";
    }
  }
  return ss.str();
}

// Electrons: local index 0 empty, 1 up, 2 dn, 3 up&dn. Up is shown as a blue
// arrow, dn as an orange arrow, double occupancy as both, empty as a dot.
static std::string to_string(ProductState const &state, Electron const &) {
  std::stringstream ss;
  for (int64_t i = state.size() - 1; i >= 0; --i) {
    bool up = state[i] & 1;
    bool dn = state[i] & 2;
    std::string s;
    if (up && dn) {
      const char *s = "\u2195";
      ss << fmt::format(fg(fmt::color::red), s);
    } else if (up) {
      const char *s = "\u2191";
      ss << fmt::format(fg(fmt::color::light_blue), s);
    } else if (dn) {
      const char *s = "\u2193";
      ss << fmt::format(fg(fmt::color::orange), s);
    } else { // empty
      const char *s = "\u25CC";
      ss << fmt::format(fg(fmt::color::gray), s);
    }
  }
  return ss.str();
}

std::string to_string(ProductState const &state, Block const &block) {
  return std::visit([&](auto const &b) { return to_string(state, b); }, block);
}

} // namespace xdiag
