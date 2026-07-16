// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::random {

uint64_t hash(Permutation const &perm);
uint64_t hash(PermutationGroup const &group);
uint64_t hash(Representation const &irrep);

uint64_t hash(Block const &block);


} // namespace xdiag::random
