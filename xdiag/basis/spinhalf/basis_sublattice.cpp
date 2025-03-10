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

template <typename bit_t, typename coeff_t, int n_sublat>
static std::pair<std::vector<bit_t>, std::vector<double>>
reps_norms_no_sz(GroupActionSublattice<bit_t, n_sublat> const &group_action,
                 arma::Col<coeff_t> const &characters) {
  using combinatorics::Subsets;

  std::vector<bit_t> reps;
  std::vector<double> norms;

  int64_t nsites = group_action.nsites();
  int64_t nsites_sublat = nsites / n_sublat;
  int64_t n_leading = nsites_sublat;
  int64_t n_trailing = (n_sublat - 1) * nsites_sublat;

  for (auto prefix : Subsets<bit_t>(n_leading)) {

    // if prefix is not rep, the full state also cannot be rep
    bit_t prefix_rep = group_action.reps_[n_sublat - 1][prefix];
    if (prefix_rep < prefix) {
      continue;
    }
#ifndef _OPENMP
    for (auto postfix : Subsets<bit_t>(n_trailing)) {
      bit_t state = (prefix << n_trailing) | postfix;
      bit_t rep = group_action.representative(state);
      if (state == rep) {
        double norm = symmetries::norm(rep, group_action, characters);
        if (std::abs(norm) > 1e-6) {
          reps.push_back(rep);
          norms.push_back(norm);
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
      for (auto postfix : combinatorics::SubsetsThread<bit_t>(n_trailing)) {
        bit_t state = (prefix << n_trailing) | postfix;
        bit_t rep = group_action.representative(state);
        if (state == rep) {
          double norm = symmetries::norm(rep, group_action, characters);
          if (std::abs(norm) > 1e-6) {
            reps_thread[myid].push_back(rep);
            norms_thread[myid].push_back(norm);
          }
        }
      }
    } // pragma omp parallel

    auto reps_prefix = omp::combine_vectors(reps_thread);
    auto norms_prefix = omp::combine_vectors(norms_thread);
    reps.insert(reps.end(), reps_prefix.begin(), reps_prefix.end());
    norms.insert(norms.end(), norms_prefix.begin(), norms_prefix.end());
#endif
  } // for (auto prefix : ...)

  reps.shrink_to_fit();
  norms.shrink_to_fit();
  return {reps, norms};
}

template <typename bit_t, typename coeff_t, int n_sublat>
static std::pair<std::vector<bit_t>, std::vector<double>>
reps_norms_sz(int64_t nup,
              GroupActionSublattice<bit_t, n_sublat> const &group_action,
              arma::Col<coeff_t> const &characters) {

  using bits::popcnt;
  std::vector<bit_t> reps;
  std::vector<double> norms;

  int64_t nsites = group_action.nsites();
  int64_t nsites_sublat = nsites / n_sublat;
  int64_t n_leading = nsites_sublat;
  int64_t n_trailing = (n_sublat - 1) * nsites_sublat;

  for (auto prefix : combinatorics::Subsets<bit_t>(n_leading)) {

    // if prefix is incompatible with nup, continue
    int64_t nup_prefix = popcnt(prefix);
    int64_t nup_postfix = nup - nup_prefix;
    if ((nup_postfix < 0) || (nup_postfix > n_trailing)) {
      continue;
    }

    // if prefix is not rep, the full state also cannot be rep
    bit_t prefix_rep = group_action.reps_[n_sublat - 1][prefix];
    if (prefix_rep < prefix) {
      continue;
    }

#ifndef _OPENMP
    for (auto postfix :
         combinatorics::Combinations<bit_t>(n_trailing, nup_postfix)) {
      bit_t state = (prefix << n_trailing) | postfix;
      bit_t rep = group_action.representative(state);
      if (state == rep) {
        double norm = symmetries::norm(rep, group_action, characters);
        if (std::abs(norm) > 1e-6) {
          reps.push_back(rep);
          norms.push_back(norm);
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
           combinatorics::CombinationsThread<bit_t>(n_trailing, nup_postfix)) {
        bit_t state = (prefix << n_trailing) | postfix;
        bit_t rep = group_action.representative(state);
        if (state == rep) {
          double norm = symmetries::norm(rep, group_action, characters);
          if (std::abs(norm) > 1e-6) {
            reps_thread[myid].push_back(rep);
            norms_thread[myid].push_back(norm);
          }
        }
      }
    } // pragma omp parallel

    auto reps_prefix = omp::combine_vectors(reps_thread);
    auto norms_prefix = omp::combine_vectors(norms_thread);
    reps.insert(reps.end(), reps_prefix.begin(), reps_prefix.end());
    norms.insert(norms.end(), norms_prefix.begin(), norms_prefix.end());
#endif

  } // for (auto prefix : ...)

  reps.shrink_to_fit();
  norms.shrink_to_fit();
  return {reps, norms};
}

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
    Representation const &irrep) try
    : nsites_(irrep.group().nsites()), nup_(undefined),
      n_postfix_bits_(nsites_ - std::min(maximum_prefix_bits, nsites_)),
      group_action_(irrep.group()), irrep_(irrep) {
  check_nsites_work_with_bits<bit_t>(nsites_);

  if (isreal(irrep)) {
    arma::vec characters = irrep.characters().as<arma::vec>();
    std::tie(reps_, norms_) = reps_norms_no_sz(group_action_, characters);
  } else {
    arma::cx_vec characters = irrep.characters().as<arma::cx_vec>();
    std::tie(reps_, norms_) = reps_norms_no_sz(group_action_, characters);
  }

  rep_search_range_ = compute_rep_search_range_serial(reps_, n_postfix_bits_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename bit_t, int n_sublat>
BasisSublattice<bit_t, n_sublat>::BasisSublattice(
    int64_t nup, Representation const &irrep) try
    : nsites_(irrep.group().nsites()), nup_(nup),
      n_postfix_bits_(nsites_ - std::min(maximum_prefix_bits, nsites_)),
      group_action_(irrep.group()), irrep_(irrep) {
  check_nsites_work_with_bits<bit_t>(nsites_);

  if (isreal(irrep)) {
    arma::vec characters = irrep.characters().as<arma::vec>();
    std::tie(reps_, norms_) = reps_norms_sz(nup, group_action_, characters);
  } else {
    arma::cx_vec characters = irrep.characters().as<arma::cx_vec>();
    std::tie(reps_, norms_) = reps_norms_sz(nup, group_action_, characters);
  }

#ifdef _OPENMP
  // omp version still has a bug, use serial (not so bad performance actually)
  // rep_search_range_ = compute_rep_search_range_omp(reps_, n_postfix_bits_);
  rep_search_range_ = compute_rep_search_range_serial(reps_, n_postfix_bits_);
#else
  rep_search_range_ = compute_rep_search_range_serial(reps_, n_postfix_bits_);
#endif
} catch (Error const &e) {
  XDIAG_RETHROW(e);
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
int64_t BasisSublattice<bit_t, n_sublat>::nsites() const {
  return nsites_;
}

template <typename bit_t, int n_sublat>
int64_t BasisSublattice<bit_t, n_sublat>::nup() const {
  return nup_;
}

template <typename bit_t, int n_sublat>

Representation const &BasisSublattice<bit_t, n_sublat>::irrep() const {
  return irrep_;
}

template <typename bit_t, int n_sublat>
GroupActionSublattice<bit_t, n_sublat> const &
BasisSublattice<bit_t, n_sublat>::group_action() const {
  return group_action_;
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
  return (nsites_ == rhs.nsites_) && (nup_ == rhs.nup_) &&
         (n_postfix_bits_ == rhs.n_postfix_bits_) &&
         (group_action_ == rhs.group_action_);
}

template <typename bit_t, int n_sublat>
bool BasisSublattice<bit_t, n_sublat>::operator!=(
    BasisSublattice<bit_t, n_sublat> const &rhs) const {
  return !operator==(rhs);
}

template class BasisSublattice<uint32_t, 1>;
template class BasisSublattice<uint32_t, 2>;
template class BasisSublattice<uint32_t, 3>;
template class BasisSublattice<uint32_t, 4>;
template class BasisSublattice<uint32_t, 5>;

template class BasisSublattice<uint64_t, 1>;
template class BasisSublattice<uint64_t, 2>;
template class BasisSublattice<uint64_t, 3>;
template class BasisSublattice<uint64_t, 4>;
template class BasisSublattice<uint64_t, 5>;

} // namespace xdiag::basis::spinhalf
