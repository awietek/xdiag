// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "mapping_syms.hpp"

namespace xdiag::symmetries {

template <typename bit_t, class group_action_t>
std::vector<int64_t> mapping_syms(bit_t origin, bit_t target,
                                  group_action_t const &group_action) {
  std::vector<int64_t> syms;
  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, origin);
    if (tstate == target) {
      syms.push_back(sym);
    }
  }
  return syms;
}

} // namespace xdiag::symmetries
