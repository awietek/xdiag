#pragma once

#include <xdiag/common.h>
#include <xdiag/symmetries/permutation.h>
#include <xdiag/symmetries/permutation_group.h>

#include <xdiag/blocks/blocks.h>
#include <xdiag/blocks/electron/electron.h>
#include <xdiag/blocks/spinhalf/spinhalf.h>
#include <xdiag/blocks/tj/tj.h>
#include <xdiag/blocks/tj_distributed/tj_distributed.h>

namespace xdiag::random {

uint64_t hash(Permutation const &perm);
uint64_t hash(PermutationGroup const &group);
uint64_t hash(Representation const &irrep);

uint64_t hash(block_variant_t const &block);
uint64_t hash(Spinhalf const &spinhalf);
uint64_t hash(tJ const &tj);
uint64_t hash(Electron const &electron);

#ifdef XDIAG_USE_MPI
uint64_t hash(tJDistributed const &tj);
#endif

} // namespace xdiag::random
