#include "block.hpp"

#include <xdiag/operators/logic/qns.hpp>
namespace xdiag {

Spinhalf block(OpSum const &ops, Spinhalf const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nupi = block.n_up();
  auto irrepi = block.irrep();
  if (!nupi && !irrepi) {
    return block;
  } else if (nupi && !irrepi) {
    auto nupr = nup(ops, block);
    return (nupi == nupr) ? block : Spinhalf(n_sites, nupr, backend);
  } else if (!nupi && irrepi) {
    auto irrepr = representation(ops, block);
    return isapprox(*irrepi, irrepr) ? block
                                     : Spinhalf(n_sites, irrepr, backend);
  } else { //(nup && irrep)
    auto nupr = nup(ops, block);
    auto irrepr = representation(ops, block);
    return ((*nupi == nupr) && isapprox(*irrepi, irrepr))
               ? block
               : Spinhalf(n_sites, nupr, irrepr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

tJ block(OpSum const &ops, tJ const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nupi = block.n_up();
  auto ndni = block.n_dn();
  auto irrepi = block.irrep();
  if (!nupi && !irrepi) {
    return block;
  } else if (nupi && !irrepi) {
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    return ((*nupi == nupr) && (*ndni == ndnr))
               ? block
               : tJ(n_sites, nupr, ndnr, backend);
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
    return ((*nupi == nupr) && (*ndni == ndnr) && isapprox(*irrepi, irrepr))
               ? block
               : tJ(n_sites, nupr, ndnr, irrepr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Electron block(OpSum const &ops, Electron const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nupi = block.n_up();
  auto ndni = block.n_dn();
  auto irrepi = block.irrep();
  if (!nupi && !irrepi) {
    return block;
  } else if (nupi && !irrepi) {
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    return ((*nupi == nupr) && (*ndni == ndnr))
               ? block
               : Electron(n_sites, nupr, ndnr, backend);
  } else if (!nupi && irrepi) {
    auto irrepr = representation(ops, block);
    return isapprox(*irrepi, irrepr) ? block
                                     : Electron(n_sites, irrepr, backend);
  } else { //(nup && irrep)
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    auto irrepr = representation(ops, block);
    return ((*nupi == nupr) && (*ndni == ndnr) && isapprox(*irrepi, irrepr))
               ? block
               : Electron(n_sites, nupr, ndnr, irrepr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

#ifdef XDIAG_USE_MPI
SpinhalfDistributed block(OpSum const &ops,
                          SpinhalfDistributed const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nupi = block.n_up();
  if (!nupi) {
    return block;
  } else {
    auto nupr = nup(ops, block);
    return (*nupi == nupr) ? block
                           : SpinhalfDistributed(n_sites, nupr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
tJDistributed block(OpSum const &ops, tJDistributed const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nupi = block.n_up();
  auto ndni = block.n_dn();
  if (!nupi) {
    return block;
  } else {
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    return ((*nupi == nupr) && (*ndni == ndnr))
               ? block
               : tJDistributed(n_sites, nupr, ndnr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
#endif

Block block(OpSum const &ops, Block const &blocki) try {
  return std::visit([&](auto &&b) { return Block(block(ops, b)); }, blocki);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool blocks_match(OpSum const &ops, Block const &block1,
                  Block const &block2) try {
  return std::visit(
      overload{
          [&](Spinhalf const &b1, Spinhalf const &b2) {
            return blocks_match(ops, b1, b2);
          },
          [&](tJ const &b1, tJ const &b2) {
	    return blocks_match(ops, b1, b2); },
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
#endif
          [&](auto const &b1, auto const &b2) { return false; },
      },
      block1, block2);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool blocks_match(OpSum const &ops, Spinhalf const &b1,
                  Spinhalf const &b2) try {
  bool match_nup = b1.n_up() ? nup(ops, b1) == *b2.n_up() : !b2.n_up();
  bool match_irrep =
      b1.irrep() ? representation(ops, b1) == *b2.irrep() : !b2.irrep();
  return match_nup && match_irrep;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool blocks_match(OpSum const &ops, tJ const &b1, tJ const &b2) try {
  bool match_nup = b1.n_up() ? nup(ops, b1) == *b2.n_up() : !b2.n_up();
  bool match_ndn = b1.n_dn() ? ndn(ops, b1) == *b2.n_dn() : !b2.n_dn();
  bool match_irrep =
      b1.irrep() ? representation(ops, b1) == *b2.irrep() : !b2.irrep();
  return match_nup && match_ndn && match_irrep;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool blocks_match(OpSum const &ops, Electron const &b1,
                  Electron const &b2) try {
  bool match_nup = b1.n_up() ? nup(ops, b1) == *b2.n_up() : !b2.n_up();
  bool match_ndn = b1.n_dn() ? ndn(ops, b1) == *b2.n_dn() : !b2.n_dn();
  bool match_irrep =
      b1.irrep() ? representation(ops, b1) == *b2.irrep() : !b2.irrep();
  return match_nup && match_ndn && match_irrep;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

#ifdef XDIAG_USE_MPI
bool blocks_match(OpSum const &ops, SpinhalfDistributed const &b1,
                  SpinhalfDistributed const &b2) try {
  bool match_nup = b1.n_up() ? nup(ops, b1) == *b2.n_up() : !b2.n_up();
  return match_nup;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool blocks_match(OpSum const &ops, tJDistributed const &b1,
                  tJDistributed const &b2) try {
  bool match_nup = b1.n_up() ? nup(ops, b1) == *b2.n_up() : !b2.n_up();
  bool match_ndn = b1.n_dn() ? ndn(ops, b1) == *b2.n_dn() : !b2.n_dn();
  return match_nup && match_ndn;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
#endif
} // namespace xdiag
