// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "random_state.hpp"

namespace xdiag {

RandomState::RandomState(int64_t seed, bool normalized)
    : seed_(seed), normalized_(normalized) {}
int64_t RandomState::seed() const { return seed_; }
bool RandomState::normalized() const { return normalized_; }

std::ostream &operator<<(std::ostream &out, RandomState const &state) {
  out << "  RandomState, seed  : " << state.seed() << "\n";
  return out;
}
std::string to_string(RandomState const &state) {
  return to_string_generic(state);
}

} // namespace xdiag
