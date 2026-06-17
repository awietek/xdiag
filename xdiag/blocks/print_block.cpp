// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

// format.hpp defines FMT_HEADER_ONLY before pulling in fmt; it must precede any
// other fmt include (print_block.hpp -> fmt/color.hpp) so the locale-aware
// format_de ("{:L}") resolves header-only instead of against an uncompiled fmt.
#include <xdiag/utils/format.hpp>

#include "print_block.hpp"

#include <cstdint>
#include <optional>
#include <regex>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/random/hash.hpp>
#include <xdiag/symmetries/representation.hpp>

#include <xdiag/extern/fmt/format.hpp>

namespace xdiag {

template <typename block_t>
void print_block(std::ostream &out, block_t const &block) {

  // Block-specific bits: name, header color, and the conserved quantum numbers.
  std::string name;
  fmt::text_style color;
  std::vector<std::pair<std::string, std::optional<int64_t>>> charges;
  if constexpr (std::is_same_v<block_t, Spinhalf>) {
    name = "Spinhalf";
    color = fmt::fg(fmt::color::moccasin) | fmt::emphasis::bold;
    charges = {{"nup", block.irreps().charge("nup")}};
  } else if constexpr (std::is_same_v<block_t, Boson>) {
    name = "Boson";
    color = fmt::fg(fmt::color::steel_blue) | fmt::emphasis::bold;
    charges = {{"number", block.irreps().charge("number")}};
  } else if constexpr (std::is_same_v<block_t, Fermion>) {
    name = "Fermion";
    color = fmt::fg(fmt::color::crimson) | fmt::emphasis::bold;
    charges = {{"number", block.irreps().charge("number")}};
  } else if constexpr (std::is_same_v<block_t, Electron>) {
    name = "Electron";
    color = fmt::fg(fmt::color::blue_violet) | fmt::emphasis::bold;
    charges = {{"nup", block.irreps().charge("nup")},
               {"ndn", block.irreps().charge("ndn")}};
  } else if constexpr (std::is_same_v<block_t, tJ>) {
    name = "tJ";
    color = fmt::fg(fmt::color::medium_sea_green) | fmt::emphasis::bold;
    charges = {{"nup", block.irreps().charge("nup")},
               {"ndn", block.irreps().charge("ndn")}};
  }

  out << fmt::format(color, "{}\n", name);
  out << fmt::format("│ {:<9}: {}\n", "nsites", block.nsites());
  out << fmt::format("│ {:<9}: {}\n", "d", block.d());
  for (auto const &[label, value] : charges) {
    if (value) {
      out << fmt::format("│ {:<9}: {}\n", label, *value);
    } else {
      out << fmt::format("│ {:<9}: not conserved\n", label);
    }
  }

  auto group = block.irreps().group("SitePermutation");
  if (group) {
    out << "│ permutation symmetries used\n";
    out << fmt::format(
        "│ {:<9}: {:x}\n", "irrep ID",
        random::hash(Representation(
            *group, *block.irreps().characters("SitePermutation"))));
  }

  std::string basisname(block.basis()->name());
  basisname = std::regex_replace(basisname, std::regex("xdiag::basis::"), "");
  basisname =
      std::regex_replace(basisname, std::regex("xdiag::combinatorics::"), "");
  basisname =
      std::regex_replace(basisname, std::regex("unsigned int"), "uint32_t");

  out << fmt::format("│ {:<9}: {}\n", "basis", basisname);
  out << fmt::format("│ {:<9}: {:x}\n", "ID", random::hash(block));
  out << fmt::format("│ {:<9}: {}\n", "dimension",
                     fmt::format_de("{:L}", block.size()));
}

template void print_block<Spinhalf>(std::ostream &, Spinhalf const &);
template void print_block<Boson>(std::ostream &, Boson const &);
template void print_block<Fermion>(std::ostream &, Fermion const &);
template void print_block<Electron>(std::ostream &, Electron const &);
template void print_block<tJ>(std::ostream &, tJ const &);

} // namespace xdiag
