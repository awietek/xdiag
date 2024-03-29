#pragma once

#include <hydra/common.h>
#include <hydra/symmetries/permutation.h>
#include <hydra/symmetries/permutation_group.h>

#include <hydra/blocks/blocks.h>
#include <hydra/blocks/electron/electron.h>
#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/tj_distributed/tj_distributed.h>

namespace hydra::random {

uint64_t hash(Permutation const &perm);
uint64_t hash(PermutationGroup const &group);
uint64_t hash(Representation const &irrep);

uint64_t hash(block_variant_t const &block);
uint64_t hash(Spinhalf const &spinhalf);
uint64_t hash(tJ const &tj);
uint64_t hash(Electron const &electron);

#ifdef HYDRA_USE_MPI
uint64_t hash(tJDistributed const &tj);
#endif

} // namespace hydra::random
