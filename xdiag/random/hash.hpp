#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/electron/electron.hpp>
#include <xdiag/blocks/spinhalf/spinhalf.hpp>
#include <xdiag/blocks/tj/tj.hpp>
#include <xdiag/blocks/tj_distributed/tj_distributed.hpp>
#include <xdiag/common.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/symmetries/permutation_group.hpp>

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
