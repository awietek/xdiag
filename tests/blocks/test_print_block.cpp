// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <sstream>
#include <string>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/utils/error.hpp>

// Exercises print_block via to_string for every block flavor. print_block has
// no coverage otherwise; here we only assert it produces a non-trivial header
// containing the fields it promises (dimension, and for symmetric blocks the
// irrep line).
TEST_CASE("print_block", "[blocks]") try {
  using namespace xdiag;

  // --- Non-symmetric blocks, with and without conserved charges ---
  {
    std::string s = to_string(Spinhalf(4));
    REQUIRE(!s.empty());
    REQUIRE(s.find("nsites") != std::string::npos);
    REQUIRE(s.find("dimension") != std::string::npos);
  }
  {
    std::string s = to_string(Spinhalf(4, 2));
    REQUIRE(s.find("dimension") != std::string::npos);
  }
  {
    std::string s = to_string(Electron(4));
    REQUIRE(s.find("dimension") != std::string::npos);
  }
  {
    std::string s = to_string(Electron(4, 2, 1));
    REQUIRE(s.find("dimension") != std::string::npos);
  }
  {
    std::string s = to_string(tJ(4, 2, 1));
    REQUIRE(s.find("dimension") != std::string::npos);
  }
  {
    std::string s = to_string(Fermion(4));
    REQUIRE(s.find("dimension") != std::string::npos);
  }
  {
    std::string s = to_string(Fermion(4, 2));
    REQUIRE(s.find("dimension") != std::string::npos);
  }
  {
    std::string s = to_string(Boson(4, 2));
    REQUIRE(s.find("dimension") != std::string::npos);
  }
  {
    std::string s = to_string(Boson(4, 2, 2));
    REQUIRE(s.find("dimension") != std::string::npos);
  }

  // --- Symmetric blocks: exercise the permutation-symmetry branch (real and
  //     complex irreps) ---
  {
    auto irrep = cyclic_group_irrep(4, 0); // trivial, real
    std::string s = to_string(Spinhalf(4, 2, irrep));
    REQUIRE(s.find("irrep ID") != std::string::npos);
    REQUIRE(s.find("dimension") != std::string::npos);
  }
  {
    auto irrep = cyclic_group_irrep(4, 1); // complex momentum
    std::string s = to_string(Spinhalf(4, 2, irrep));
    REQUIRE(s.find("irrep ID") != std::string::npos);
  }
  {
    auto irrep = cyclic_group_irrep(4, 0);
    std::string s = to_string(Electron(4, 2, 1, irrep));
    REQUIRE(s.find("irrep ID") != std::string::npos);
  }
  {
    auto irrep = cyclic_group_irrep(4, 0);
    std::string s = to_string(tJ(4, 2, 1, irrep));
    REQUIRE(s.find("irrep ID") != std::string::npos);
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("block product-state printing", "[blocks]") try {
  using namespace xdiag;

  // to_string(ProductState, Block) is the per-block pretty-printer. Exercise it
  // for every flavor with a valid product state, plus operator<<(Block).
  {
    auto block = Spinhalf(4);
    ProductState p(std::vector<int64_t>{1, 0, 1, 0});
    std::string s = to_string(p, Block(block));
    REQUIRE(!s.empty());

    std::ostringstream oss;
    oss << Block(block);
    REQUIRE(!oss.str().empty());
  }
  {
    auto block = Electron(4);
    ProductState p(std::vector<int64_t>{0, 1, 2, 3});
    REQUIRE(!to_string(p, Block(block)).empty());
  }
  {
    auto block = tJ(4, 2, 1);
    ProductState p(std::vector<int64_t>{1, 2, 0, 1});
    REQUIRE(!to_string(p, Block(block)).empty());
  }
  {
    auto block = Boson(4, 2);
    ProductState p(std::vector<int64_t>{0, 1, 1, 0});
    REQUIRE(!to_string(p, Block(block)).empty());
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
