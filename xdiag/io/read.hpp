#pragma once

#include <string>

#include <xdiag/common.hpp>

#include <xdiag/io/file_toml.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

XDIAG_API Permutation read_permutation(FileToml file, std::string tag);
XDIAG_API PermutationGroup read_permutation_group(FileToml file,
                                                  std::string tag);
XDIAG_API Representation read_representation(
    FileToml file, std::string irrep_tag, std::string group_tag = "Symmetries");
XDIAG_API Op read_op(FileToml file, std::string tag);
XDIAG_API OpSum read_opsum(FileToml file, std::string tag);

} // namespace xdiag
