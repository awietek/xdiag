// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "representative.hpp"

#include <xdiag/bits/bitset.hpp>

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

template std::pair<uint32_t, std::vector<int64_t>>
representative_syms<uint32_t>(uint32_t, SitePermutation const &);
template std::pair<uint64_t, std::vector<int64_t>>
representative_syms<uint64_t>(uint64_t, SitePermutation const &);

template <typename bit_t>
bit_t representative_subset(bit_t state, SitePermutation const &action,
                            std::vector<int64_t> const &syms) {
  bit_t rep = state;
  bool first = true;
  for (int64_t sym : syms) {
    bit_t trans = action.apply(sym, state);
    if (first || (trans < rep)) {
      rep = trans;
      first = false;
    }
  }
  return rep;
}

template uint32_t representative_subset<uint32_t>(uint32_t,
                                                  SitePermutation const &,
                                                  std::vector<int64_t> const &);
template uint64_t representative_subset<uint64_t>(uint64_t,
                                                  SitePermutation const &,
                                                  std::vector<int64_t> const &);
template bits::BitsetDynamic representative_subset<bits::BitsetDynamic>(
    bits::BitsetDynamic, SitePermutation const &, std::vector<int64_t> const &);

template <typename bit_t>
std::pair<bit_t, int64_t>
representative_sym_subset(bit_t state, SitePermutation const &action,
                          std::vector<int64_t> const &syms) {
  bit_t rep = state;
  int64_t rep_sym = 0;
  bool first = true;
  for (int64_t sym : syms) {
    bit_t trans = action.apply(sym, state);
    if (first || (trans < rep)) {
      rep = trans;
      rep_sym = sym;
      first = false;
    }
  }
  return {rep, rep_sym};
}

template std::pair<uint32_t, int64_t>
representative_sym_subset<uint32_t>(uint32_t, SitePermutation const &,
                                    std::vector<int64_t> const &);
template std::pair<uint64_t, int64_t>
representative_sym_subset<uint64_t>(uint64_t, SitePermutation const &,
                                    std::vector<int64_t> const &);
template std::pair<bits::BitsetDynamic, int64_t>
representative_sym_subset<bits::BitsetDynamic>(bits::BitsetDynamic,
                                               SitePermutation const &,
                                               std::vector<int64_t> const &);

} // namespace xdiag::symmetries
