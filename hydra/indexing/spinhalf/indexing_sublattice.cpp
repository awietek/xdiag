#include "indexing_sublattice.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>

#include <hydra/symmetries/operations/group_action_operations.h>
#include <hydra/symmetries/operations/symmetry_operations.h>
#include <hydra/utils/logger.h>

#include <algorithm>

#ifdef _OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#endif

namespace hydra::indexing::spinhalf {

template <typename bit_t>
ska::flat_hash_map<bit_t, gsl::span<bit_t const>>
compute_rep_search_range_serial(std::vector<bit_t> const &reps,
                                int n_postfix_bits) {
  ska::flat_hash_map<bit_t, gsl::span<bit_t const>> rep_search_range;

  if (reps.size() > 0) {
    idx_t start = 0;
    idx_t end = 0;
    bit_t previous_prefix = reps[0] >> n_postfix_bits;
    for (auto rep : reps) {
      bit_t prefix = rep >> n_postfix_bits;
      if (prefix != previous_prefix) {
        rep_search_range[previous_prefix] =
            gsl::span<bit_t const>(reps.data() + start, end - start);
        previous_prefix = prefix;
        start = end;
      }
      ++end;
    }
    rep_search_range[previous_prefix] =
        gsl::span<bit_t const>(reps.data() + start, end - start);
  }

  return rep_search_range;
}

#ifdef _OPENMP
template <typename bit_t>
ska::flat_hash_map<bit_t, gsl::span<bit_t const>>
compute_rep_search_range_omp(std::vector<bit_t> const &reps,
                             int n_postfix_bits) {
  //// COMMENT: HAS A BUG DO NOT USE

  ska::flat_hash_map<bit_t, gsl::span<bit_t const>> rep_search_range;

  if (reps.size() > 0) {

    std::vector<std::vector<std::pair<bit_t, gsl::span<bit_t const>>>>
        rep_search_range_thread;

#pragma omp parallel
    {
      int myid = omp_get_thread_num();
      int rank = omp_get_num_threads();

#pragma omp single
      { rep_search_range_thread.resize(rank); }

      auto [start_thread, end_thread] = omp::get_omp_start_end(reps.size());

      // Adjust start/end such that there is no overlaps between threads
      if ((myid != rank - 1) && (end_thread != 0) &&
          (end_thread != (idx_t)reps.size())) {
        while (reps[end_thread - 1] == reps[end_thread]) {
          ++end_thread;
        }
      }

      if ((myid != 0) && (start_thread != 0) &&
          (end_thread != (idx_t)reps.size())) {
        while (reps[start_thread + 1] == reps[start_thread]) {
          ++start_thread;
        }
      }

      idx_t start = start_thread;
      idx_t end = start_thread;
      bit_t previous_prefix = reps[start] >> n_postfix_bits;
      for (idx_t idx = start_thread; idx < end_thread; ++idx) {
        bit_t prefix = reps[idx] >> n_postfix_bits;
        if (prefix != previous_prefix) {

          rep_search_range_thread[myid].push_back(
              {previous_prefix,
               gsl::span<bit_t const>(reps.data() + start, end - start)});

          previous_prefix = prefix;
          start = end;
        }
        ++end;
      }
      rep_search_range_thread[myid].push_back(
          {previous_prefix,
           gsl::span<bit_t const>(reps.data() + start, end - start)});
    } // #pragma omp parallel

    for (auto const &prefixes_ranges : rep_search_range_thread) {
      for (auto [prefix, range] : prefixes_ranges) {
        rep_search_range[prefix] = range;
      }
    }
  }

  return rep_search_range;
} // namespace hydra::indexing::spinhalf
#endif

template <typename bit_t, int n_sublat>
IndexingSublattice<bit_t, n_sublat>::IndexingSublattice(
    int n_sites, PermutationGroup permutation_group, Representation irrep)
    : n_sites_(n_sites), sz_conserved_(false), n_up_(invalid_n),
      n_postfix_bits_(n_sites - std::min(maximum_prefix_bits, n_sites)),
      group_action_(allowed_subgroup(permutation_group, irrep)), irrep_(irrep) {
  using combinatorics::Subsets;

  int n_sites_sublat = n_sites_ / n_sublat;
  int n_leading = n_sites_sublat;
  int n_trailing = (n_sublat - 1) * n_sites_sublat;

  for (auto prefix : Subsets<bit_t>(n_leading)) {

    // if prefix is not rep, the full state also cannot be rep
    bit_t prefix_rep = group_action_.reps_[n_sublat - 1][prefix];
    if (prefix_rep < prefix) {
      continue;
    }

    for (auto postfix : Subsets<bit_t>(n_trailing)) {
      bit_t state = (prefix << n_trailing) | postfix;
      bit_t rep = representative(state);
      if (state == rep) {
        double norm = symmetries::norm(rep, group_action_, irrep);
        if (std::abs(norm) > 1e-6) {
          reps_.push_back(rep);
          norms_.push_back(norm);
        }
      }
    }
  }

  reps_.shrink_to_fit();
  norms_.shrink_to_fit();
  rep_search_range_ = compute_rep_search_range_serial(reps_, n_postfix_bits_);
  begin_ = iterator_t(reps_, 0);
  end_ = iterator_t(reps_, reps_.size());
}

template <typename bit_t, int n_sublat>
IndexingSublattice<bit_t, n_sublat>::IndexingSublattice(
    int n_sites, int n_up, PermutationGroup permutation_group,
    Representation irrep)
    : n_sites_(n_sites), sz_conserved_(true), n_up_(n_up),
      n_postfix_bits_(n_sites - std::min(maximum_prefix_bits, n_sites)),
      group_action_(allowed_subgroup(permutation_group, irrep)), irrep_(irrep) {
  using bitops::popcnt;

  int n_sites_sublat = n_sites_ / n_sublat;
  int n_leading = n_sites_sublat;
  int n_trailing = (n_sublat - 1) * n_sites_sublat;

  for (auto prefix : combinatorics::Subsets<bit_t>(n_leading)) {

    // if prefix is incompatible with n_up, continue
    int n_up_prefix = popcnt(prefix);
    int n_up_postfix = n_up - n_up_prefix;
    if ((n_up_postfix < 0) || (n_up_postfix > n_trailing)) {
      continue;
    }

    // if prefix is not rep, the full state also cannot be rep
    bit_t prefix_rep = group_action_.reps_[n_sublat - 1][prefix];
    if (prefix_rep < prefix) {
      continue;
    }

#ifndef _OPENMP
    for (auto postfix :
         combinatorics::Combinations<bit_t>(n_trailing, n_up_postfix)) {
      bit_t state = (prefix << n_trailing) | postfix;
      bit_t rep = representative(state);
      if (state == rep) {
        double norm = symmetries::norm(rep, group_action_, irrep);
        if (std::abs(norm) > 1e-6) {
          reps_.push_back(rep);
          norms_.push_back(norm);
        }
      }
    }

#else

    std::vector<std::vector<bit_t>> reps_thread;
    std::vector<std::vector<double>> norms_thread;
#pragma omp parallel
    {
      int myid = omp_get_thread_num();
      int rank = omp_get_num_threads();
#pragma omp single
      {
        reps_thread.resize(rank);
        norms_thread.resize(rank);
      }

      // Compute representatives for each thread
      for (auto postfix :
           combinatorics::CombinationsThread<bit_t>(n_trailing, n_up_postfix)) {
        bit_t state = (prefix << n_trailing) | postfix;
        bit_t rep = representative(state);
        if (state == rep) {
          double norm = symmetries::norm(rep, group_action_, irrep);
          if (std::abs(norm) > 1e-6) {
            reps_thread[myid].push_back(rep);
            norms_thread[myid].push_back(norm);
          }
        }
      }
    } // pragma omp parallel

    auto reps_prefix = omp::combine_vectors(reps_thread);
    auto norms_prefix = omp::combine_vectors(norms_thread);
    reps_.insert(reps_.end(), reps_prefix.begin(), reps_prefix.end());
    norms_.insert(norms_.end(), norms_prefix.begin(), norms_prefix.end());
#endif

  } // for (auto prefix : ...)

  reps_.shrink_to_fit();
  norms_.shrink_to_fit();
#ifdef _OPENMP
  // omp version still has a bug, use serial (not so bad performance actually)
  // rep_search_range_ = compute_rep_search_range_omp(reps_, n_postfix_bits_);
  rep_search_range_ = compute_rep_search_range_serial(reps_, n_postfix_bits_);
#else
  rep_search_range_ = compute_rep_search_range_serial(reps_, n_postfix_bits_);
#endif

  begin_ = iterator_t(reps_, 0);
  end_ = iterator_t(reps_, reps_.size());
} // namespace hydra::indexing::spinhalf

template <typename bit_t, int n_sublat>
idx_t IndexingSublattice<bit_t, n_sublat>::index_of_representative(
    bit_t rep) const {

  bit_t prefix = rep >> n_postfix_bits_;
  auto itr = rep_search_range_.find(prefix);
  if (itr == rep_search_range_.end()) {
    return invalid_index;
  } else {
    gsl::span<bit_t const> search_range = itr->second;
    auto it = std::lower_bound(search_range.begin(), search_range.end(), rep);
    if ((it != search_range.end()) && (*it == rep)) {
      return &(*it) - reps_.data();
    } else {
      return invalid_index;
    }
  }
}

template <typename bit_t, int n_sublat>
idx_t IndexingSublattice<bit_t, n_sublat>::index(bit_t state) const {
  bit_t rep = representative(state);
  return index_of_representative(rep);
}

template <typename bit_t, int n_sublat>
std::pair<idx_t, int>
IndexingSublattice<bit_t, n_sublat>::index_sym(bit_t state) const {
  auto [rep, sym] = group_action_.representative_sym(state);
  idx_t idx = index_of_representative(rep);
  return {idx, sym};
}

template <typename bit_t, int n_sublat>
std::pair<idx_t, gsl::span<int const>>
IndexingSublattice<bit_t, n_sublat>::index_syms(bit_t state) const {
  auto [rep, syms] = group_action_.representative_syms(state);
  idx_t idx = index_of_representative(rep);
  return {idx, syms};
}

template <typename bit_t, int n_sublat>
bool IndexingSublattice<bit_t, n_sublat>::operator==(
    IndexingSublattice<bit_t, n_sublat> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (n_postfix_bits_ == rhs.n_postfix_bits_) &&
         (group_action_ == rhs.group_action_) && (irrep_ == rhs.irrep_);
}

template <typename bit_t, int n_sublat>
bool IndexingSublattice<bit_t, n_sublat>::operator!=(
    IndexingSublattice<bit_t, n_sublat> const &rhs) const {
  return !operator==(rhs);
}

template class IndexingSublattice<uint16_t, 1>;
template class IndexingSublattice<uint32_t, 1>;
template class IndexingSublattice<uint64_t, 1>;

template class IndexingSublattice<uint16_t, 2>;
template class IndexingSublattice<uint32_t, 2>;
template class IndexingSublattice<uint64_t, 2>;

template class IndexingSublattice<uint16_t, 3>;
template class IndexingSublattice<uint32_t, 3>;
template class IndexingSublattice<uint64_t, 3>;

template class IndexingSublattice<uint16_t, 4>;
template class IndexingSublattice<uint32_t, 4>;
template class IndexingSublattice<uint64_t, 4>;

template class IndexingSublattice<uint16_t, 5>;
template class IndexingSublattice<uint32_t, 5>;
template class IndexingSublattice<uint64_t, 5>;

} // namespace hydra::indexing::spinhalf
