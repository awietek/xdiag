#include "read.hpp"

#include <xdiag/utils/vector.hpp>

namespace xdiag {
Permutation read_permutation(FileToml file, std::string tag) try {
  if (!file.defined(tag)) {
    XDIAG_THROW(fmt::format("TOML file does not contain the tag \"{}\"", tag));
  }
  return file[tag].as<Permutation>();
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

PermutationGroup read_permutation_group(FileToml file, std::string tag) try {
  if (!file.defined(tag)) {
    XDIAG_THROW(fmt::format("TOML file does not contain the tag \"{}\"", tag));
  }
  return file[tag].as<PermutationGroup>();
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
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
    return Representation(sgroup, characters);
  } else {
    return Representation(group, characters);
  }

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Op read_op(FileToml file, std::string tag) try {
  if (!file.defined(tag)) {
    XDIAG_THROW(fmt::format("TOML file does not contain the tag \"{}\"", tag));
  }
  return file[tag].as<Op>();
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
OpSum read_opsum(FileToml file, std::string tag) try {
  if (!file.defined(tag)) {
    XDIAG_THROW(fmt::format("TOML file does not contain the tag \"{}\"", tag));
  }
  return file[tag].as<OpSum>();
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
} // namespace xdiag
