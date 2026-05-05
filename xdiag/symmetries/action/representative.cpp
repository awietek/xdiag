// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "representative.hpp"

namespace xdiag::symmetries {

template <typename bit_t>
bit_t representative(bit_t state, SitePermutation const &action) {
  bit_t rep = std::numeric_limits<bit_t>::max();
  for (int64_t sym = 0; sym < action.size(); ++sym) {
    bit_t trans = action.apply(sym, state);
    if (trans < rep) {
      rep = trans;
    }
  }
  return rep;
}

template uint16_t representative<uint16_t>(uint16_t, SitePermutation const &);
template uint32_t representative<uint32_t>(uint32_t, SitePermutation const &);
template uint64_t representative<uint64_t>(uint64_t, SitePermutation const &);

template <typename bit_t>
std::pair<bit_t, int64_t> representative_sym(bit_t state,
                                             SitePermutation const &action) {
  bit_t rep = std::numeric_limits<bit_t>::max();
  int64_t idx = 0;
  for (int64_t sym = 0; sym < action.size(); ++sym) {
    bit_t trans = action.apply(sym, state);
    if (trans < rep) {
      rep = trans;
      idx = sym;
    }
  }
  return {rep, idx};
}

template std::pair<uint16_t, int64_t>
representative_sym<uint16_t>(uint16_t, SitePermutation const &);
template std::pair<uint32_t, int64_t>
representative_sym<uint32_t>(uint32_t, SitePermutation const &);
template std::pair<uint64_t, int64_t>
representative_sym<uint64_t>(uint64_t, SitePermutation const &);

template <typename bit_t>
std::pair<bit_t, std::vector<int64_t>>
representative_syms(bit_t state, SitePermutation const &action) {
  bit_t rep = std::numeric_limits<bit_t>::max();
  std::vector<int64_t> indices;
  for (int64_t sym = 0; sym < action.size(); ++sym) {
    bit_t trans = action.apply(sym, state);
    if (trans < rep) {
      rep = trans;
      indices.clear();
      indices.push_back(sym);
    } else if (trans == rep) {
      indices.push_back(sym);
    }
  }
  return {rep, indices};
}

template std::pair<uint16_t, std::vector<int64_t>>
representative_syms<uint16_t>(uint16_t, SitePermutation const &);
template std::pair<uint32_t, std::vector<int64_t>>
representative_syms<uint32_t>(uint32_t, SitePermutation const &);
template std::pair<uint64_t, std::vector<int64_t>>
representative_syms<uint64_t>(uint64_t, SitePermutation const &);

} // namespace xdiag::symmetries
