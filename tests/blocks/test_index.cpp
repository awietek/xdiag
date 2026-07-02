// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/utils/error.hpp>

using namespace xdiag;

// For every basis state produced by iterating a block, both the member
// block.index(pstate) and the free index(Block, pstate) must return that
// state's iteration position -- i.e. index is the exact O(1) inverse of the
// ProductState iterator. Also checks the count matches size().
template <typename block_t>
static void check_roundtrip(block_t const &block) {
  int64_t i = 0;
  for (ProductState const &pstate : block) {
    REQUIRE(block.index(pstate) == i);         // concrete-block member
    REQUIRE(index(Block(block), pstate) == i); // free fn over the Block variant
    ++i;
  }
  REQUIRE(i == block.size());
}

TEST_CASE("index_spinhalf", "[blocks]") try {
  for (int64_t N = 1; N <= 8; ++N) {
    check_roundtrip(Spinhalf(N)); // full Hilbert space
    for (int64_t nup = 0; nup <= N; ++nup) {
      check_roundtrip(Spinhalf(N, nup));
    }
  }
  // translational symmetry, every momentum sector
  for (int64_t N : {2, 4, 6}) {
    for (int64_t k = 0; k < N; ++k) {
      Representation irrep = cyclic_group_irrep(N, k);
      check_roundtrip(Spinhalf(N, irrep));
      for (int64_t nup = 0; nup <= N; ++nup) {
        check_roundtrip(Spinhalf(N, nup, irrep));
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("index_fermion", "[blocks]") try {
  for (int64_t N = 1; N <= 8; ++N) {
    check_roundtrip(Fermion(N));
    for (int64_t nf = 0; nf <= N; ++nf) {
      check_roundtrip(Fermion(N, nf));
    }
  }
  for (int64_t N : {2, 4, 6}) {
    for (int64_t k = 0; k < N; ++k) {
      Representation irrep = cyclic_group_irrep(N, k);
      check_roundtrip(Fermion(N, irrep));
      for (int64_t nf = 0; nf <= N; ++nf) {
        check_roundtrip(Fermion(N, nf, irrep));
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("index_electron", "[blocks]") try {
  for (int64_t N = 1; N <= 5; ++N) {
    check_roundtrip(Electron(N));
    for (int64_t nup = 0; nup <= N; ++nup) {
      for (int64_t ndn = 0; ndn <= N; ++ndn) {
        check_roundtrip(Electron(N, nup, ndn));
      }
    }
  }
  for (int64_t N : {2, 4}) {
    for (int64_t k = 0; k < N; ++k) {
      Representation irrep = cyclic_group_irrep(N, k);
      check_roundtrip(Electron(N, irrep));
      for (int64_t nup = 0; nup <= N; ++nup) {
        for (int64_t ndn = 0; ndn <= N; ++ndn) {
          check_roundtrip(Electron(N, nup, ndn, irrep));
        }
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("index_tj", "[blocks]") try {
  // t-J: no double occupancy, so nup + ndn <= N
  for (int64_t N = 1; N <= 5; ++N) {
    for (int64_t nup = 0; nup <= N; ++nup) {
      for (int64_t ndn = 0; nup + ndn <= N; ++ndn) {
        check_roundtrip(tJ(N, nup, ndn));
      }
    }
  }
  for (int64_t N : {2, 4}) {
    for (int64_t k = 0; k < N; ++k) {
      Representation irrep = cyclic_group_irrep(N, k);
      for (int64_t nup = 0; nup <= N; ++nup) {
        for (int64_t ndn = 0; nup + ndn <= N; ++ndn) {
          check_roundtrip(tJ(N, nup, ndn, irrep));
        }
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("index_boson", "[blocks]") try {
  // d = local dimension (max occupation = d - 1)
  for (int64_t N = 1; N <= 4; ++N) {
    for (int64_t d = 2; d <= 4; ++d) {
      check_roundtrip(Boson(N, d)); // number not conserved
      for (int64_t nb = 0; nb <= (d - 1) * N; ++nb) {
        check_roundtrip(Boson(N, d, nb));
      }
    }
  }
  for (int64_t N : {2, 3, 4}) {
    for (int64_t k = 0; k < N; ++k) {
      Representation irrep = cyclic_group_irrep(N, k);
      for (int64_t d = 2; d <= 3; ++d) {
        check_roundtrip(Boson(N, d, irrep));
        for (int64_t nb = 0; nb <= (d - 1) * N; ++nb) {
          check_roundtrip(Boson(N, d, nb, irrep));
        }
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}
