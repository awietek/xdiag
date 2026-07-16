// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <ostream>
#include <string>

#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

class XDIAG_API RandomState {
public:
  RandomState(int64_t seed = 42, bool normalized = true);
  int64_t seed() const;
  bool normalized() const;

private:
  int64_t seed_;
  bool normalized_;
};

XDIAG_API std::ostream &operator<<(std::ostream &out, RandomState const &state);
XDIAG_API std::string to_string(RandomState const &state);

} // namespace xdiag
