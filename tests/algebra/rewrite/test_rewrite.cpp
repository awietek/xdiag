// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <optional>
#include <vector>

#include <tests/catch.hpp>

#include <xdiag/algebra/rewrite/rewrite.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace xdiag::algebra;

// The engine is a single generic fixed-point driver: rewrite_to_fixpoint applies
// a one-step rewrite to every monomial, collects, and repeats until stable.
TEST_CASE("test_rewrite", "[operators]") try {
  Log("Testing rewrite_to_fixpoint");

  // A step that lowers an "N" operator to "Cdag" "C", and is a no-op otherwise.
  // Applied to a fixed point, an OpSum of N's becomes Cdag*C products.
  auto lower_n = [](Monomial const &mono) -> std::optional<OpSum> {
    for (int64_t k = 0; k < mono.size(); ++k) {
      if (mono[k].type() != "N") {
        continue;
      }
      int64_t i = mono[k][0];
      std::vector<Op> repl;
      for (int64_t j = 0; j < k; ++j) {
        repl.push_back(mono[j]);
      }
      repl.push_back(Op("Cdag", i));
      repl.push_back(Op("C", i));
      for (int64_t j = k + 1; j < mono.size(); ++j) {
        repl.push_back(mono[j]);
      }
      return OpSum(Monomial(repl));
    }
    return std::nullopt;
  };

  {
    // 2.0 * N{0}  ->  2.0 * Cdag{0} C{0}
    OpSum r = rewrite_to_fixpoint(2.0 * Op("N", 0), lower_n);
    REQUIRE(r == 2.0 * (Op("Cdag", 0) * Op("C", 0)));
  }
  {
    // N{0} * N{1}  ->  Cdag{0} C{0} Cdag{1} C{1}  (both expanded)
    OpSum r = rewrite_to_fixpoint(1.0 * (Op("N", 0) * Op("N", 1)), lower_n);
    REQUIRE(r ==
            1.0 * (Op("Cdag", 0) * Op("C", 0) * Op("Cdag", 1) * Op("C", 1)));
  }
  {
    // A step returning an empty OpSum drops the monomial (collected to zero).
    auto annihilate = [](Monomial const &mono) -> std::optional<OpSum> {
      if (mono.size() >= 1 && mono[0].type() == "Z") {
        return OpSum{};
      }
      return std::nullopt;
    };
    OpSum r = rewrite_to_fixpoint(3.0 * Op("Z", 0), annihilate);
    REQUIRE(r == OpSum{});
  }

  Log("Done testing rewrite_to_fixpoint");
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
