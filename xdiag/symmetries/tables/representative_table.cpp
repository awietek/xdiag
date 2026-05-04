// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "representative_table.hpp"

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <set>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/log2.hpp>
#include <xdiag/bits/nbits.hpp>
#include <xdiag/combinatorics/bounded_multisets/bounded_multisets.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>
#include <xdiag/combinatorics/bounded_partitions/schaefer_table.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/symmetries/action/isrepresentative.hpp>
#include <xdiag/symmetries/action/norm.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::symmetries {

template <typename enumeration_t, typename coeff_t>
static void representative_table_initialize(
    enumeration_t const &enumeration, SitePermutation const &action,
    arma::Col<coeff_t> const &characters,
    bits::BitVector<typename enumeration_t::bit_t> &representative,
    bits::BitVector<uint64_t> &representative_index,
    bits::BitVector<uint64_t> &representative_symmetry,
    bits::BitVector<uint64_t> &representative_norm_index,
    std::vector<double> &norms) try {
  using bit_t = typename enumeration_t::bit_t;
  using bits::BitVector;

  if (action.size() != characters.size()) {
    XDIAG_THROW(fmt::format("Size of symmetry group action (={}) incompatible "
                            "with size of character table of irrep (={})",
                            action.size(), characters.size()));
  }
  if (action.nsites() != enumeration.n()) {
    XDIAG_THROW(
        fmt::format("nsites of symmetry group action (={}) incompatible "
                    "with nsites of state enumeration (={})",
                    action.nsites(), enumeration.n()));
  }

  // First pass, simply count number of representatives and number of different
  // norms
#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
  std::vector<int64_t> nrepresentatives_for_thread(nthreads, 0);
  std::vector<std::set<double>> norms_for_thread(nthreads);

#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    int64_t size = enumeration.size();
    auto begin = enumeration.begin() + num_thread * (size / nthreads);
    auto end = (num_thread == nthreads - 1)
                   ? enumeration.end()
                   : enumeration.begin() + (num_thread + 1) * (size / nthreads);
    for (auto it = begin; it != end; ++it) {
      auto state = *it;
      if (isrepresentative(state, action)) {
        double nrm = norm(state, action, characters);
        if (std::fabs(nrm) > 1e-6) { // representative found
          ++nrepresentatives_for_thread[num_thread];
          norms_for_thread[num_thread].insert(nrm);
        }
      }
    }
  }

  // combine norms from different threads
  std::set<double> norms_set;
  for (auto const &s : norms_for_thread) {
    norms_set.insert(s.begin(), s.end());
  }
  norms = std::vector<double>(norms_set.begin(), norms_set.end());

  // get the offsets for nrepresentatives
  std::vector<int64_t> nrepresentatives_for_thread_offset(nthreads, 0);
  std::partial_sum(nrepresentatives_for_thread.begin(),
                   nrepresentatives_for_thread.end() - 1,
                   nrepresentatives_for_thread_offset.begin() + 1);

  // get the total number of representatives
  int64_t nrepresentatives =
      std::accumulate(nrepresentatives_for_thread.begin(),
                      nrepresentatives_for_thread.end(), 0);
#else
  int64_t nrepresentatives = 0;
  for (auto state : enumeration) {
    if (isrepresentative(state, action)) {
      double nrm = norm(state, action, characters);

      if (std::fabs(nrm) > 1e-6) { // representative found
        ++nrepresentatives;

        // Check if norm already registered ...
        auto it = std::find_if(norms.begin(), norms.end(), [&](double n) {
          return std::fabs(n - nrm) < 1e-6;
        });

        // ... if not, register it
        if (it == norms.end()) {
          norms.push_back(nrm);
        }
      }
    }
  }
#endif
  std::sort(norms.begin(), norms.end());

  // --------------------------------------------------------------
  // Now come the allocations, since we know al the relevant sizes
  // --------------------------------------------------------------

  // Create vector holding the representatives
  try {
    int64_t size = nrepresentatives;
    int64_t nbits = enumeration.bitwidth();
    representative = BitVector<bit_t>(size, nbits);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative array");
  }

  // Create vector holding the indices for each state yielding the
  // representative
  try {
    int64_t size = enumeration.size();
    int64_t nbits = std::max(1u, bits::ceillog2(nrepresentatives + 1));
    representative_index = BitVector<uint64_t>(size, nbits);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative index array");
  }

  // Create vector holding the symmetry which yields the representative
  try {
    int64_t size = enumeration.size();
    int64_t nbits = std::max(1u, bits::ceillog2(action.size()));
    representative_symmetry = BitVector<uint64_t>(size, nbits);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative symmetry array");
  }

  // Create vector holding the norm index of the states
  try {
    int64_t size = nrepresentatives;
    int64_t nbits = std::max(1u, bits::ceillog2(norms.size()));
    representative_norm_index = BitVector<uint64_t>(size, nbits);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative norm index array");
  }

  // --------------------------------------------------------------
  // Second pass: fill representative[] and representative_norm_index[].
  // Thread t writes to the contiguous range [rep_offset[t], rep_offset[t+1]).
  // Adjacent thread ranges may share a 64-bit storage chunk in the BitVector,
  // so we use atomic_or_element (OR into zero-initialised storage) which is
  // safe because each element position is written by exactly one thread.
  // --------------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel
  {
    int t = omp_get_thread_num();
    int64_t size = enumeration.size();
    int64_t i_begin = (int64_t)t * (size / nthreads);
    int64_t i_end =
        (t == nthreads - 1) ? size : (int64_t)(t + 1) * (size / nthreads);
    auto it = enumeration.begin() + i_begin;
    int64_t local_rep_idx = nrepresentatives_for_thread_offset[t];
    for (int64_t i = i_begin; i < i_end; ++i, ++it) {
      bit_t state = *it;
      if (isrepresentative(state, action)) {
        double nrm = norm(state, action, characters);
        if (std::fabs(nrm) > 1e-6) {
          representative.atomic_or_element(local_rep_idx, state);
          auto it2 = std::find_if(norms.begin(), norms.end(), [&](double n) {
            return std::fabs(n - nrm) < 1e-6;
          });
          representative_norm_index.atomic_or_element(
              local_rep_idx, (uint64_t)std::distance(norms.begin(), it2));
          ++local_rep_idx;
        }
      }
    }
  }
#else
  {
    int64_t rep_idx = 0;
    for (auto state : enumeration) {
      if (isrepresentative(state, action)) {
        double nrm = norm(state, action, characters);
        if (std::fabs(nrm) > 1e-6) {
          representative[rep_idx] = state;
          auto it = std::find_if(norms.begin(), norms.end(), [&](double n) {
            return std::fabs(n - nrm) < 1e-6;
          });
          representative_norm_index[rep_idx] =
              (uint64_t)std::distance(norms.begin(), it);
          ++rep_idx;
        }
      }
    }
  }
#endif

  // --------------------------------------------------------------
  // Third pass: expand each representative's orbit to fill
  // representative_index[] and representative_symmetry[] for every state.
  // Orbits are disjoint so each idx is written by exactly one rep_idx.
  // When a representative has a non-trivial stabilizer, multiple symmetries
  // map it to the same orbit member; we deduplicate via a per-thread `seen`
  // list so each idx receives exactly one write (first symmetry wins).
  // In the OMP path we use atomic_or_element (OR into zero-initialised
  // storage) so concurrent writes to different rep_idx are race-free.
  // --------------------------------------------------------------
  auto const &group = action.group();
#ifdef _OPENMP
#pragma omp parallel
  {
    std::vector<int64_t> seen;
    seen.reserve(action.size());
#pragma omp for schedule(static)
    for (int64_t rep_idx = 0; rep_idx < (int64_t)representative.size();
         ++rep_idx) {
      seen.clear();
      bit_t rep = representative[rep_idx];
      for (int64_t sym = 0; sym < action.size(); ++sym) {
        bit_t state = action.apply(sym, rep);
        int64_t idx = enumeration.index(state);
        if (std::find(seen.begin(), seen.end(), idx) != seen.end())
          continue;
        seen.push_back(idx);
        representative_index.atomic_or_element(idx, (uint64_t)(rep_idx + 1));
        representative_symmetry.atomic_or_element(idx,
                                                  (uint64_t)group.inv(sym));
      }
    }
  }
#else
  for (int64_t rep_idx = 0; rep_idx < (int64_t)representative.size();
       ++rep_idx) {
    bit_t rep = representative[rep_idx];
    for (int64_t sym = 0; sym < action.size(); ++sym) {
      bit_t state = action.apply(sym, rep);
      int64_t idx = enumeration.index(state);
      representative_index[idx] = (uint64_t)(rep_idx + 1);
      representative_symmetry[idx] = (uint64_t)group.inv(sym);
    }
  }
#endif
}
XDIAG_CATCH

template <typename enumeration_t>
RepresentativeTable<enumeration_t>::RepresentativeTable(
    enumeration_t const &enumeration, Representation const &irrep) try {
  auto action = SitePermutation(irrep.group());
  if (isreal(irrep)) {
    representative_table_initialize(
        enumeration, action, irrep.characters().as<arma::vec>(),
        representative_, representative_index_, representative_symmetry_,
        representative_norm_index_, norm_);
  } else {
    representative_table_initialize(
        enumeration, action, irrep.characters().as<arma::cx_vec>(),
        representative_, representative_index_, representative_symmetry_,
        representative_norm_index_, norm_);
  }
}
XDIAG_CATCH

template <typename enumeration_t>
int64_t RepresentativeTable<enumeration_t>::size() const {
  return representative_.size();
}

template <typename enumeration_t>
bool RepresentativeTable<enumeration_t>::operator==(
    RepresentativeTable const &rhs) const {
  return (representative_ == rhs.representative_) &&
         (representative_index_ == rhs.representative_index_) &&
         (representative_symmetry_ == rhs.representative_symmetry_) &&
         (representative_norm_index_ == rhs.representative_norm_index_) &&
         (norm_ == rhs.norm_);
}

template <typename enumeration_t>
bool RepresentativeTable<enumeration_t>::operator!=(
    RepresentativeTable const &rhs) const {
  return !operator==(rhs);
}

} // namespace xdiag::symmetries

using namespace xdiag::combinatorics;
using namespace xdiag::bits;

#define INSTANTIATE_REPRESENTATIVE_TABLE(ENUMERATION_TYPE)                     \
  template class xdiag::symmetries::RepresentativeTable<ENUMERATION_TYPE>;

// BEGIN_INSTANTIATION_GROUP(subsets)
INSTANTIATE_REPRESENTATIVE_TABLE(Subsets<uint32_t>);
INSTANTIATE_REPRESENTATIVE_TABLE(Subsets<uint64_t>);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(combinations)
INSTANTIATE_REPRESENTATIVE_TABLE(Combinations<uint32_t>);
INSTANTIATE_REPRESENTATIVE_TABLE(Combinations<uint64_t>);
INSTANTIATE_REPRESENTATIVE_TABLE(Combinations<BitsetStatic2>);
INSTANTIATE_REPRESENTATIVE_TABLE(Combinations<BitsetStatic4>);
INSTANTIATE_REPRESENTATIVE_TABLE(Combinations<BitsetStatic8>);
INSTANTIATE_REPRESENTATIVE_TABLE(Combinations<BitsetDynamic>);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(lintable)
INSTANTIATE_REPRESENTATIVE_TABLE(LinTable<uint32_t>);
INSTANTIATE_REPRESENTATIVE_TABLE(LinTable<uint64_t>);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(bounded_multisets)
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray1>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray2>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray3>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray4>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray5>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray6>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray7>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray8>);

INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong1>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong2>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong3>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong4>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong5>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong6>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong7>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong8>);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(bounded_partitions)
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray1>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray2>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray3>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray4>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray5>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray6>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray7>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray8>);

INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong1>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong2>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong3>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong4>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong5>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong6>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong7>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong8>);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(schaefer_table)
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray1>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray2>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray3>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray4>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray5>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray6>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray7>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray8>);
// END_INSTANTIATION_GROUP

#undef INSTANTIATE_REPRESENTATIVE_TABLE
