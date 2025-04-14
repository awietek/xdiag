// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/common.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

class RandomState {
public:
  XDIAG_API RandomState(int64_t seed = 42, bool normalized = true);
  int64_t seed() const;
  bool normalized() const;

private:
  int64_t seed_;
  bool normalized_;
};

XDIAG_API std::ostream &operator<<(std::ostream &out, RandomState const &state);
XDIAG_API std::string to_string(RandomState const &state);

} // namespace xdiag
