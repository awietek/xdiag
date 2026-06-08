// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/utils/logger.hpp>

TEST_CASE("boson", "[boson]") try {
  using namespace xdiag;
  int64_t nsites = 4;
  int64_t d = 6;
  auto block = Boson(nsites, d, d*2);
  for (auto state : block) {
    Log(to_string(state, block));
  }
} catch (xdiag::Error e) {
  error_trace(e);
}
