// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "block.hpp"

#include <xdiag/operators/logic/order.hpp>
#include <xdiag/operators/logic/qns.hpp>
#include <xdiag/utils/xdiag_show.hpp>

namespace xdiag {

Spinhalf block(OpSum const &ops, Spinhalf const &block) try {
  int64_t nsites = block.nsites();
  std::string backend = block.backend();
  auto nupi = block.nup();
  auto irrepi = block.irrep();
  if (!nupi && !irrepi) {
    return block;
  } else if (nupi && !irrepi) {
    auto nupr = nup(ops, block);
    return (nupi == nupr) ? block : Spinhalf(nsites, nupr, backend);
  } else if (!nupi && irrepi) {
    auto irrepr = representation(ops, block);
    return isapprox(*irrepi, irrepr) ? block
                                     : Spinhalf(nsites, irrepr, backend);
  } else { //(nup && irrep)
    auto nupr = nup(ops, block);
    auto irrepr = representation(ops, block);
    return ((*nupi == nupr) && isapprox(*irrepi, irrepr))
               ? block
               : Spinhalf(nsites, nupr, irrepr, backend);
  }
}
XDIAG_CATCH

tJ block(OpSum const &ops, tJ const &block) try {
  int64_t nsites = block.nsites();
  std::string backend = block.backend();
  auto nupi = block.nup();
  auto ndni = block.ndn();
  auto irrepi = block.irrep();
  if (!nupi && !irrepi) {
    return block;
  } else if (nupi && !irrepi) {
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    return ((*nupi == nupr) && (*ndni == ndnr))
               ? block
               : tJ(nsites, nupr, ndnr, backend);
  }
  // // Not yet implemented
  // else if (!nup && irrep) {
  //   auto irrepr = representation(ops, block);
  //   return isapprox(irrep, irrepr) ? block : tJ(nsites, *irrep, backend);
  // }
  else { //(nup && irrep)
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    auto irrepr = representation(ops, block);
    return ((*nupi == nupr) && (*ndni == ndnr) && isapprox(*irrepi, irrepr))
               ? block
               : tJ(nsites, nupr, ndnr, irrepr, backend);
  }
}
XDIAG_CATCH

Electron block(OpSum const &ops, Electron const &block) try {
  int64_t nsites = block.nsites();
  std::string backend = block.backend();
  auto nupi = block.nup();
  auto ndni = block.ndn();
  auto irrepi = block.irrep();
  if (!nupi && !irrepi) {
    return block;
  } else if (nupi && !irrepi) {
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    return ((*nupi == nupr) && (*ndni == ndnr))
               ? block
               : Electron(nsites, nupr, ndnr, backend);
  } else if (!nupi && irrepi) {
    auto irrepr = representation(ops, block);
    return isapprox(*irrepi, irrepr) ? block
                                     : Electron(nsites, irrepr, backend);
  } else { //(nup && irrep)
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    auto irrepr = representation(ops, block);
    return ((*nupi == nupr) && (*ndni == ndnr) && isapprox(*irrepi, irrepr))
               ? block
               : Electron(nsites, nupr, ndnr, irrepr, backend);
  }
}
XDIAG_CATCH

#ifdef XDIAG_USE_MPI
SpinhalfDistributed block(OpSum const &ops,
                          SpinhalfDistributed const &block) try {
  int64_t nsites = block.nsites();
  std::string backend = block.backend();
  auto nupi = block.nup();
  if (!nupi) {
    return block;
  } else {
    auto nupr = nup(ops, block);
    return (*nupi == nupr) ? block : SpinhalfDistributed(nsites, nupr, backend);
  }
}
XDIAG_CATCH
tJDistributed block(OpSum const &ops, tJDistributed const &block) try {
  int64_t nsites = block.nsites();
  std::string backend = block.backend();
  auto nupi = block.nup();
  auto ndni = block.ndn();
  if (!nupi) {
    return block;
  } else {
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    return ((*nupi == nupr) && (*ndni == ndnr))
               ? block
               : tJDistributed(nsites, nupr, ndnr, backend);
  }
}
XDIAG_CATCH

ElectronDistributed block(OpSum const &ops,
                          ElectronDistributed const &block) try {
  int64_t nsites = block.nsites();
  std::string backend = block.backend();
  auto nupi = block.nup();
  auto ndni = block.ndn();
  if (!nupi) {
    return block;
  } else {
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    return ((*nupi == nupr) && (*ndni == ndnr))
               ? block
               : ElectronDistributed(nsites, nupr, ndnr, backend);
  }
}
XDIAG_CATCH
#endif

Block block(OpSum const &ops, Block const &blocki) try {
  return std::visit([&](auto &&b) { return Block(block(ops, b)); }, blocki);
}
XDIAG_CATCH

bool blocks_match(OpSum const &ops, Block const &block1,
                  Block const &block2) try {
  return std::visit(
      overload{
          [&](Spinhalf const &b1, Spinhalf const &b2) {
            return blocks_match(ops, b1, b2);
          },
          [&](tJ const &b1, tJ const &b2) { return blocks_match(ops, b1, b2); },
          [&](Electron const &b1, Electron const &b2) {
            return blocks_match(ops, b1, b2);
          },
#ifdef XDIAG_USE_MPI
          [&](SpinhalfDistributed const &b1, SpinhalfDistributed const &b2) {
            return blocks_match(ops, b1, b2);
          },
          [&](tJDistributed const &b1, tJDistributed const &b2) {
            return blocks_match(ops, b1, b2);
          },
          [&](ElectronDistributed const &b1, ElectronDistributed const &b2) {
            return blocks_match(ops, b1, b2);
          },
#endif
          [&](auto const &b1, auto const &b2) { return false; },
      },
      block1, block2);
}
XDIAG_CATCH

bool blocks_match(OpSum const &ops, Spinhalf const &b1,
                  Spinhalf const &b2) try {
  bool match_nup = b1.nup() ? nup(ops, b1) == *b2.nup() : !b2.nup();
  bool match_irrep =
      b1.irrep() ? isapprox(representation(ops, b1), *b2.irrep()) : !b2.irrep();
  return match_nup && match_irrep;
}
XDIAG_CATCH

bool blocks_match(OpSum const &ops, tJ const &b1, tJ const &b2) try {
  bool match_nup = b1.nup() ? nup(ops, b1) == *b2.nup() : !b2.nup();
  bool match_ndn = b1.ndn() ? ndn(ops, b1) == *b2.ndn() : !b2.ndn();
  bool match_irrep =
      b1.irrep() ? isapprox(representation(ops, b1), *b2.irrep()) : !b2.irrep();
  return match_nup && match_ndn && match_irrep;
}
XDIAG_CATCH

bool blocks_match(OpSum const &ops, Electron const &b1,
                  Electron const &b2) try {
  bool match_nup = b1.nup() ? nup(ops, b1) == *b2.nup() : !b2.nup();
  bool match_ndn = b1.ndn() ? ndn(ops, b1) == *b2.ndn() : !b2.ndn();
  bool match_irrep =
      b1.irrep() ? isapprox(representation(ops, b1), *b2.irrep()) : !b2.irrep();
  return match_nup && match_ndn && match_irrep;
}
XDIAG_CATCH

#ifdef XDIAG_USE_MPI
bool blocks_match(OpSum const &ops, SpinhalfDistributed const &b1,
                  SpinhalfDistributed const &b2) try {
  bool match_nup = b1.nup() ? nup(ops, b1) == *b2.nup() : !b2.nup();
  return match_nup;
}
XDIAG_CATCH

bool blocks_match(OpSum const &ops, tJDistributed const &b1,
                  tJDistributed const &b2) try {
  bool match_nup = b1.nup() ? nup(ops, b1) == *b2.nup() : !b2.nup();
  bool match_ndn = b1.ndn() ? ndn(ops, b1) == *b2.ndn() : !b2.ndn();
  return match_nup && match_ndn;
}
XDIAG_CATCH

bool blocks_match(OpSum const &ops, ElectronDistributed const &b1,
                  ElectronDistributed const &b2) try {
  bool match_nup = b1.nup() ? nup(ops, b1) == *b2.nup() : !b2.nup();
  bool match_ndn = b1.ndn() ? ndn(ops, b1) == *b2.ndn() : !b2.ndn();
  return match_nup && match_ndn;
}
XDIAG_CATCH
#endif
} // namespace xdiag
