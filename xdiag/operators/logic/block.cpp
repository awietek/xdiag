#include "block.hpp"

namespace xdiag {

Spinhalf block(OpSum const &ops, Spinhalf const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nup = block.nup();
  auto irrep = block.irrep();
  if (!nup && !irrep) {
    return block;
  } else if (nup && !irrep) {
    auto nupr = nup(ops, block);
    return (nup == nupr) ? block : Spinhalf(n_sites, *nupr, backend);
  } else if (!nup && irrep) {
    auto irrepr = representation(ops, block);
    return isapprox(irrep, irrepr) ? block : Spinhalf(n_sites, *irrep, backend);
  } else { //(nup && irrep)
    auto nupr = nup(ops, block);
    auto irrepr = representation(ops, block);
    return ((nup == nupr) && isapprox(irrep, irrepr))
               ? block
               : Spinhalf(n_sites, *nupr, *irrepr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

tJ block(OpSum const &ops, tJ const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nup = block.nup();
  auto ndn = block.ndn();
  auto irrep = block.irrep();
  if (!nup && !irrep) {
    return block;
  } else if (nup && !irrep) {
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    return ((nup == nupr) && (ndn == ndnr))
               ? block
               : tJ(n_sites, *nupr, *ndnr, backend);
  }
  // // Not yet implemented
  // else if (!nup && irrep) {
  //   auto irrepr = representation(ops, block);
  //   return isapprox(irrep, irrepr) ? block : tJ(n_sites, *irrep, backend);
  // }
  else { //(nup && irrep)
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    auto irrepr = representation(ops, block);
    return ((nup == nupr) && (ndn == ndnr) && isapprox(irrep, irrepr))
               ? block
               : tJ(n_sites, *nupr, *ndnr, *irrepr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Electron block(OpSum const &ops, Electron const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nup = block.nup();
  auto ndn = block.ndn();
  auto irrep = block.irrep();
  if (!nup && !irrep) {
    return block;
  } else if (nup && !irrep) {
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    return ((nup == nupr) && (ndn == ndnr))
               ? block
               : Electron(n_sites, *nupr, *ndnr, backend);
  } else if (!nup && irrep) {
    auto irrepr = representation(ops, block);
    return isapprox(irrep, irrepr) ? block : Electron(n_sites, *irrep, backend);
  } else { //(nup && irrep)
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    auto irrepr = representation(ops, block);
    return ((nup == nupr) && (ndn == ndnr) && isapprox(irrep, irrepr))
               ? block
               : Electron(n_sites, *nupr, *ndnr, *irrepr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

#ifdef XDIAG_USE_MPI
SpinhalfDistributed block(OpSum const &ops,
                          SpinhalfDistributed const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nup = block.nup();
  if (!nup) {
    return block;
  } else {
    auto nupr = nup(ops, block);
    return (nup == nupr) ? block : SpinhalfDistributed(n_sites, *nupr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
tJDistributed block(OpSum const &ops, tJDistributed const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nup = block.nup();
  auto ndn = block.ndn();
  if (!nup) {
    return block;
  } else {
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    return ((nup == nupr) && (ndn == ndnr))
               ? block
               : tJDistributed(n_sites, *nupr, *ndnr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
#endif

Block block(OpSum const &ops, Block const &block) try {
  return std::visit([&](auto &&b) { return block(ops, b); }, block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool blocks_match(OpSum const &ops, Block const &block1,
                  Block const &block2) try {
  return std::visit(
      overload{
          [&](Spinhalf const &b1, Spinhalf const &b2) {
            return blocks_match(b1, b2);
          },
          [&](tJ const &b1, tJ const &b2) { return blocks_match(b1, b2); },
          [&](Electron const &b1, Electron const &b2) {
            return blocks_match(b1, b2);
          },
#ifdef XDIAG_USE_MPI
          [&](SpinhalfDistributed const &b1, SpinhalfDistributed const &b2) {
            return blocks_match(b1, b2);
          },
          [&](tJDistributed const &b1, tJDistributed const &b2) {
            return blocks_match(b1, b2);
          },
#endif
          [&](auto const &b1, auto const &b2) { return false; },
      },
      block1, block2);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool blocks_match(OpSum const &ops, Spinhalf const &b1,
                  Spinhalf const &b2) try {
  auto nupr = nup(ops, b1);
  auto nup = b2.nup();
  auto irrepr = representation(ops, b1);
  auto irrep = b2.irrep();
  return (nup == nupr) && (irrep == irrepr);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool blocks_match(OpSum const &ops, tJ const &b1, tJ const &b2) try {
  auto nupr = nup(ops, b1);
  auto nup = b2.nup();
  auto ndnr = nup(ops, b1);
  auto ndn = b2.nup();
  auto irrepr = representation(ops, b1);
  auto irrep = b2.irrep();
  return (nup == nupr) && (ndn == ndnr) && (irrep == irrepr);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool blocks_match(OpSum const &ops, Electron const &b1,
                  Electron const &b2) try {
  auto nupr = nup(ops, b1);
  auto nup = b2.nup();
  auto ndnr = nup(ops, b1);
  auto ndn = b2.nup();
  auto irrepr = representation(ops, b1);
  auto irrep = b2.irrep();
  return (nup == nupr) && (ndn == ndnr) && (irrep == irrepr);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

#ifdef XDIAG_USE_MPI
bool blocks_match(OpSum const &ops, SpinhalfDistributed const &b1,
                  SpinhalfDistributed const &b2) try {
  auto nupr = nup(ops, b1);
  auto nup = b2.nup();
  return (nup == nupr);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool blocks_match(OpSum const &ops, tJDistributed const &b1,
                  tJDistributed const &b2) try {
  auto nupr = nup(ops, b1);
  auto nup = b2.nup();
  auto ndnr = nup(ops, b1);
  auto ndn = b2.nup();
  return (nup == nupr) && (ndn == ndnr);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
#endif
