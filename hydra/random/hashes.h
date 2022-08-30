#pragma once

#include <hydra/common.h>
#include <hydra/symmetries/permutation.h>
#include <hydra/symmetries/permutation_group.h>

#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/electron/electron.h>

namespace hydra::random {

uint32_t hash(Permutation const &perm);
uint32_t hash(PermutationGroup const &group);
uint32_t hash(Representation const &irrep);

template <typename bit_t> uint32_t hash(Spinhalf<bit_t> const &spinhalf);
template <typename bit_t> uint32_t hash(tJ<bit_t> const &tj);
template <typename bit_t> uint32_t hash(Electron<bit_t> const &electron);

} // namespace hydra::random
