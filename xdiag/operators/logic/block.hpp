#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

XDIAG_API Block block(OpSum const &ops, Block const &block);
XDIAG_API Spinhalf block(OpSum const &ops, Spinhalf const &block);
XDIAG_API tJ block(OpSum const &ops, tJ const &block);
XDIAG_API Electron block(OpSum const &ops, Electron const &block);
#ifdef XDIAG_USE_MPI
XDIAG_API SpinhalfDistributed block(OpSum const &ops,
                                    SpinhalfDistributed const &block);
XDIAG_API tJDistributed block(OpSum const &ops, tJDistributed const &block);
#endif

bool blocks_match(OpSum const &ops, Block const &b1, Block const &b2);
bool blocks_match(OpSum const &ops, Spinhalf const &b1, Spinhalf const &b2);
bool block(OpSum const &ops, tJ const &b1, tJ const &b2);
bool blocks_match(OpSum const &ops, Electron const &b1, Electron const &b2);
#ifdef XDIAG_USE_MPI
bool blocks_match(OpSum const &ops, SpinhalfDistributed const &b1,
                  SpinhalfDistributed const &b2);
bool blocks_match(OpSum const &ops, tJDistributed const &b1,
                  tJDistributed const &b2);
#endif

} // namespace xdiag
