// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "mapping_syms.hpp"

namespace xdiag::symmetries {

template <typename bit_t, class group_action_t>
std::vector<int64_t> stabilizer_symmetries(bit_t bits,
                                           group_action_t const &group_action) {
  std::vector<int64_t> stable_syms;
  for (int64_t sym = 0; sym < group_action.nsymmetries(); ++sym) {
    if (group.apply(sym, bits) == bits) {
      stable_syms.push_back(sym);
    }
  }
  return stable_syms;
}
} // namespace xdiag::symmetries
