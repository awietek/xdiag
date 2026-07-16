// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fermi_table.hpp"

#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/symmetries/fermi_sign.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::symmetries {

template <typename enumeration_t>
FermiTable<enumeration_t>::FermiTable(enumeration_t const &enumeration,
                                      PermutationGroup const &group) try
    : enumeration_(enumeration), size_(enumeration.size()) {
  if (enumeration.n() != group.nsites()) {
    XDIAG_THROW("nsites of the enumeration does not match the nsites of the "
                "PermutationGroup");
  }
  int64_t nsyms = group.size();
  table_.resize(nsyms * size_);

  // Filled serially: std::vector<bool> is bit-packed, so concurrent writes to
  // distinct entries that share a storage word would race (the reason the
  // earlier OpenMP version had to build per-thread vectors and concatenate). The
  // build is one-time and O(|group| * size * nsites) with the O(nsites)
  // fermi_bool_of_permutation, dwarfed by the matvecs the O(1) lookup then
  // speeds up. If it ever needs threads, give each symmetry a word-aligned row
  // (pad size to a multiple of 64) and parallelise over symmetries.
  for (int64_t sym = 0; sym < nsyms; ++sym) {
    Permutation perm = group[sym];
    int64_t idx = 0;
    for (bit_t state : enumeration_) {
      table_[sym * size_ + idx] = fermi_bool_of_permutation(state, perm);
      ++idx;
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::symmetries

using namespace xdiag;
using namespace xdiag::combinatorics;
using namespace xdiag::bits;

#define INSTANTIATE_FERMI_TABLE(ENUMERATION)                                   \
  template class xdiag::symmetries::FermiTable<ENUMERATION>;

INSTANTIATE_FERMI_TABLE(Subsets<uint32_t>)
INSTANTIATE_FERMI_TABLE(Subsets<uint64_t>)
INSTANTIATE_FERMI_TABLE(LinTable<uint32_t>)
INSTANTIATE_FERMI_TABLE(LinTable<uint64_t>)
INSTANTIATE_FERMI_TABLE(Combinations<uint32_t>)
INSTANTIATE_FERMI_TABLE(Combinations<uint64_t>)
INSTANTIATE_FERMI_TABLE(Combinations<BitsetDynamic>)

#undef INSTANTIATE_FERMI_TABLE
