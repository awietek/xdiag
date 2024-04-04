#include "basis_sublattice.hpp"

#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/subsets.hpp>

#include <xdiag/symmetries/operations/group_action_operations.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/utils/logger.hpp>

#include <algorithm>

#ifdef _OPENMP
#include <xdiag/parallel/omp/omp_utils.hpp>
#endif

namespace xdiag::basis::spinhalf {

template <typename bit_t>
ska::flat_hash_map<bit_t, gsl::span<bit_t const>>
compute_rep_search_range_serial(std::vector<bit_t> const &reps,
                                int n_postfix_bits) {
  ska::flat_hash_map<bit_t, gsl::span<bit_t const>> rep_search_range;

  if (reps.size() > 0) {
    int64_t start = 0;
    int64_t end = 0;
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
                             int64_t n_postfix_bits) {
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
          (end_thread != (int64_t)reps.size())) {
        while (reps[end_thread - 1] == reps[end_thread]) {
          ++end_thread;
        }
      }

      if ((myid != 0) && (start_thread != 0) &&
          (end_thread != (int64_t)reps.size())) {
        while (reps[start_thread + 1] == reps[start_thread]) {
          ++start_thread;
        }
      }

      int64_t start = start_thread;
      int64_t end = start_thread;
      bit_t previous_prefix = reps[start] >> n_postfix_bits;
      for (int64_t idx = start_thread; idx < end_thread; ++idx) {
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
} // namespace xdiag::basis::spinhalf
#endif

template <typename bit_t, int n_sublat>
BasisSublattice<bit_t, n_sublat>::BasisSublattice(
    int64_t n_sites, PermutationGroup permutation_group, Representation irrep)
    : n_sites_(n_sites), sz_conserved_(false), n_up_(invalid_n),
      n_postfix_bits_(n_sites - std::min(maximum_prefix_bits, n_sites)),
      group_action_(allowed_subgroup(permutation_group, irrep)), irrep_(irrep) {
  using combinatorics::Subsets;

  int64_t n_sites_sublat = n_sites_ / n_sublat;
  int64_t n_leading = n_sites_sublat;
  int64_t n_trailing = (n_sublat - 1) * n_sites_sublat;

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
}

template <typename bit_t, int n_sublat>
BasisSublattice<bit_t, n_sublat>::BasisSublattice(
    int64_t n_sites, int64_t n_up, PermutationGroup permutation_group,
    Representation irrep)
    : n_sites_(n_sites), sz_conserved_(true), n_up_(n_up),
      n_postfix_bits_(n_sites - std::min(maximum_prefix_bits, n_sites)),
      group_action_(allowed_subgroup(permutation_group, irrep)), irrep_(irrep) {
  using bits::popcnt;

  int64_t n_sites_sublat = n_sites_ / n_sublat;
  int64_t n_leading = n_sites_sublat;
  int64_t n_trailing = (n_sublat - 1) * n_sites_sublat;

  for (auto prefix : combinatorics::Subsets<bit_t>(n_leading)) {

    // if prefix is incompatible with n_up, continue
    int64_t n_up_prefix = popcnt(prefix);
    int64_t n_up_postfix = n_up - n_up_prefix;
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
}

template <typename bit_t, int n_sublat>
typename std::vector<bit_t>::const_iterator
BasisSublattice<bit_t, n_sublat>::begin() const {
  return reps_.begin();
}

template <typename bit_t, int n_sublat>
typename std::vector<bit_t>::const_iterator
BasisSublattice<bit_t, n_sublat>::end() const {
  return reps_.end();
}

template <typename bit_t, int n_sublat>
int64_t BasisSublattice<bit_t, n_sublat>::dim() const {
  return size();
}

template <typename bit_t, int n_sublat>
int64_t BasisSublattice<bit_t, n_sublat>::size() const {
  return reps_.size();
}

template <typename bit_t, int n_sublat>
int64_t BasisSublattice<bit_t, n_sublat>::n_sites() const {
  return n_sites_;
}

template <typename bit_t, int n_sublat>
bool BasisSublattice<bit_t, n_sublat>::sz_conserved() const {
  return sz_conserved_;
}

template <typename bit_t, int n_sublat>
int64_t BasisSublattice<bit_t, n_sublat>::n_up() const {
  return n_up_;
}

template <typename bit_t, int n_sublat>
GroupActionSublattice<bit_t, n_sublat> const &
BasisSublattice<bit_t, n_sublat>::group_action() const {
  return group_action_;
}

template <typename bit_t, int n_sublat>
Representation const &BasisSublattice<bit_t, n_sublat>::irrep() const {
  return irrep_;
}

template <typename bit_t, int n_sublat>
int64_t
BasisSublattice<bit_t, n_sublat>::index_of_representative(bit_t rep) const {

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
int64_t BasisSublattice<bit_t, n_sublat>::index(bit_t state) const {
  bit_t rep = representative(state);
  return index_of_representative(rep);
}

template <typename bit_t, int n_sublat>
std::pair<int64_t, int64_t>
BasisSublattice<bit_t, n_sublat>::index_sym(bit_t state) const {
  auto [rep, sym] = group_action_.representative_sym(state);
  int64_t idx = index_of_representative(rep);
  return {idx, sym};
}

template <typename bit_t, int n_sublat>
std::pair<int64_t, gsl::span<int64_t const>>
BasisSublattice<bit_t, n_sublat>::index_syms(bit_t state) const {
  auto [rep, syms] = group_action_.representative_syms(state);
  int64_t idx = index_of_representative(rep);
  return {idx, syms};
}

template <typename bit_t, int n_sublat>
bool BasisSublattice<bit_t, n_sublat>::operator==(
    BasisSublattice<bit_t, n_sublat> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (n_postfix_bits_ == rhs.n_postfix_bits_) &&
         (group_action_ == rhs.group_action_) && (irrep_ == rhs.irrep_);
}

template <typename bit_t, int n_sublat>
bool BasisSublattice<bit_t, n_sublat>::operator!=(
    BasisSublattice<bit_t, n_sublat> const &rhs) const {
  return !operator==(rhs);
}

template class BasisSublattice<uint16_t, 1>;
template class BasisSublattice<uint32_t, 1>;
template class BasisSublattice<uint64_t, 1>;

template class BasisSublattice<uint16_t, 2>;
template class BasisSublattice<uint32_t, 2>;
template class BasisSublattice<uint64_t, 2>;

template class BasisSublattice<uint16_t, 3>;
template class BasisSublattice<uint32_t, 3>;
template class BasisSublattice<uint64_t, 3>;

template class BasisSublattice<uint16_t, 4>;
template class BasisSublattice<uint32_t, 4>;
template class BasisSublattice<uint64_t, 4>;

template class BasisSublattice<uint16_t, 5>;
template class BasisSublattice<uint32_t, 5>;
template class BasisSublattice<uint64_t, 5>;

} // namespace xdiag::basis::spinhalf
