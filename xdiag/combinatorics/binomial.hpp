// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>

namespace xdiag::combinatorics {

int64_t binomial(int64_t n, int64_t k);

constexpr int64_t binom(int64_t n, int64_t k) {
  if (k > n || k < 0)
    return 0;
  int64_t res = 1;
  for (int64_t i = 1; i <= k; i++)
    res = (res * (n - i + 1)) / i;
  return res;
}

} // namespace xdiag::combinatorics
