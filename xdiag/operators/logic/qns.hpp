#pragma once

#include <optional>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

XDIAG_API Spinhalf block(OpSum const &ops, Spinhalf const &block);
XDIAG_API tJ block(OpSum const &ops, tJ const &block);
XDIAG_API Electron block(OpSum const &ops, Electron const &block);
#ifdef XDIAG_USE_MPI
XDIAG_API SpinhalfDistributed block(OpSum const &ops,
                                    SpinhalfDistributed const &block);
XDIAG_API tJDistributed block(OpSum const &ops, tJDistributed const &block);
#endif
XDIAG_API Block block(OpSum const &ops, Block const &block);

template <typename block_t>
XDIAG_API Representation representation(OpSum const &ops, block_t const &block);

XDIAG_API Representation representation(OpSum const &ops, Block const &block);
XDIAG_API Representation representation(OpSum const &ops, State const &v);

template <typename block_t>
XDIAG_API int64_t nup(OpSum const &ops, block_t const &block);

XDIAG_API int64_t nup(OpSum const &ops, Block const &block);
XDIAG_API int64_t nup(OpSum const &ops, State const &v);

template <typename block_t>
XDIAG_API int64_t ndn(OpSum const &ops, block_t const &block);

XDIAG_API int64_t ndn(OpSum const &ops, Block const &block);
XDIAG_API int64_t ndn(OpSum const &ops, State const &v);

std::optional<int64_t> nup(Op const &op);
std::optional<int64_t> nup(OpSum const &ops);
std::optional<int64_t> ndn(Op const &op);
std::optional<int64_t> ndn(OpSum const &ops);
std::optional<Representation> representation(OpSum const &ops,
                                             PermutationGroup const &group);

} // namespace xdiag
