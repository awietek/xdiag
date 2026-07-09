// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "read.hpp"

#include <xdiag/math/vector.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag {
Permutation read_permutation(FileToml file, std::string tag) try {
  if (!file.defined(tag)) {
    XDIAG_THROW(fmt::format("TOML file does not contain the tag \"{}\"", tag));
  }
  return file[tag].as<Permutation>();
}
XDIAG_CATCH

PermutationGroup read_permutation_group(FileToml file, std::string tag) try {
  if (!file.defined(tag)) {
    XDIAG_THROW(fmt::format("TOML file does not contain the tag \"{}\"", tag));
  }
  return file[tag].as<PermutationGroup>();
}
XDIAG_CATCH

Representation read_representation(FileToml file, std::string irrep_tag,
                                   std::string group_tag) try {
  if (!file.defined(group_tag)) {
    XDIAG_THROW(
        fmt::format("TOML file does not contain the tag \"{}\"", group_tag));
  }
  auto group = read_permutation_group(file, group_tag);
  std::string character_tag = irrep_tag + std::string(".characters");
  if (!file.defined(group_tag)) {
    XDIAG_THROW(fmt::format("TOML file does not contain the tag \"{}\"",
                            character_tag));
  }
  if (!file.defined(character_tag)) {
    XDIAG_THROW(fmt::format("TOML file does not contain the tag \"{}\"",
                            character_tag));
  }

  Vector characters;
  try {
    characters = Vector(file[character_tag].as<arma::vec>());
  } catch (Error const &e) {
    try {
      characters = Vector(file[character_tag].as<arma::cx_vec>());
    } catch (Error const &e) {
      XDIAG_THROW("Unable to parse characters");
    }
  }

  std::string allowed_syms_tag = irrep_tag + std::string(".allowed_symmetries");
  if (file.defined(allowed_syms_tag)) {
    auto allowed_syms = file[allowed_syms_tag].as<std::vector<int64_t>>();
    if (allowed_syms.size() != characters.size()) {
      XDIAG_THROW("Allowed symmetries and characters do not have the same "
                  "number of elements");
    }
    auto sgroup = subgroup(group, allowed_syms);
    if (isreal(characters)) {
      return Representation(sgroup, characters.as<arma::vec>());
    } else {
      return Representation(sgroup, characters.as<arma::cx_vec>());
    }
  } else {
    if (isreal(characters)) {
      return Representation(group, characters.as<arma::vec>());
    } else {
      return Representation(group, characters.as<arma::cx_vec>());
    }
  }
}
XDIAG_CATCH

Op read_op(FileToml file, std::string tag) try {
  if (!file.defined(tag)) {
    XDIAG_THROW(fmt::format("TOML file does not contain the tag \"{}\"", tag));
  }
  return file[tag].as<Op>();
}
XDIAG_CATCH

OpSum read_opsum(FileToml file, std::string tag) try {
  if (!file.defined(tag)) {
    XDIAG_THROW(fmt::format("TOML file does not contain the tag \"{}\"", tag));
  }
  return file[tag].as<OpSum>();
}
XDIAG_CATCH

} // namespace xdiag
