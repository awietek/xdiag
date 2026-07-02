// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <string>

#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
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
}
