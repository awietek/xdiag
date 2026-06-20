// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_sublattice.hpp"

#include <algorithm>
#include <limits>

#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/symmetries/action/norm.hpp>
#include <xdiag/symmetries/action/site_permutation_sublattice.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag::basis {

namespace {

constexpr int64_t undefined = std::numeric_limits<int64_t>::min();

template <typename bit_t, typename coeff_t, int n_sublat>
std::pair<std::vector<bit_t>, std::vector<double>> reps_norms_no_sz(
    symmetries::SitePermutationSublattice<bit_t, n_sublat> const &action,
    arma::Col<coeff_t> const &characters) {
  using combinatorics::Subsets;

  std::vector<bit_t> reps;
  std::vector<double> norms;

  int64_t nsites = action.nsites();
  int64_t nsites_sublat = nsites / n_sublat;
  int64_t n_leading = nsites_sublat;
  int64_t n_trailing = (n_sublat - 1) * nsites_sublat;

  for (auto prefix : Subsets<bit_t>(n_leading)) {

    // if prefix is not rep, the full state also cannot be rep
    auto prefix_rep = action.reps_[n_sublat - 1][(int64_t)prefix];
    if (prefix_rep < prefix) {
      continue;
    }

#ifndef _OPENMP
    for (auto postfix : Subsets<bit_t>(n_trailing)) {
      bit_t state = (prefix << n_trailing) | postfix;
      bit_t rep = action.representative(state);
      if (state == rep) {
        double norm = symmetries::norm(rep, action, characters);
        if (std::abs(norm) > 1e-6) {
          reps.push_back(rep);
          norms.push_back(norm);
        }
      }
    }
#else
    std::vector<std::vector<bit_t>> reps_thread;
    std::vector<std::vector<double>> norms_thread;
    auto postfixes = Subsets<bit_t>(n_trailing);
    int64_t sz = postfixes.size();

#pragma omp parallel
    {
      int myid = omp_get_thread_num();
      int nthr = omp_get_num_threads();
#pragma omp single
      {
        reps_thread.resize(nthr);
        norms_thread.resize(nthr);
      }
      int64_t start = myid * (sz / nthr);
      int64_t end = (myid == nthr - 1) ? sz : (myid + 1) * (sz / nthr);
      for (auto it = postfixes.begin() + start, eit = postfixes.begin() + end;
           it != eit; ++it) {
        bit_t postfix = *it;
        bit_t state = (prefix << n_trailing) | postfix;
        bit_t rep = action.representative(state);
        if (state == rep) {
          double norm = symmetries::norm(rep, action, characters);
          if (std::abs(norm) > 1e-6) {
            reps_thread[myid].push_back(rep);
            norms_thread[myid].push_back(norm);
          }
        }
      }
    } // omp parallel

    for (auto const &rv : reps_thread)
      reps.insert(reps.end(), rv.begin(), rv.end());
    for (auto const &nv : norms_thread)
      norms.insert(norms.end(), nv.begin(), nv.end());
#endif
  } // for prefix

  reps.shrink_to_fit();
  norms.shrink_to_fit();
  return {reps, norms};
}

template <typename bit_t, typename coeff_t, int n_sublat>
std::pair<std::vector<bit_t>, std::vector<double>> reps_norms_sz(
    int64_t nup,
    symmetries::SitePermutationSublattice<bit_t, n_sublat> const &action,
    arma::Col<coeff_t> const &characters) {
  std::vector<bit_t> reps;
  std::vector<double> norms;

  int64_t nsites = action.nsites();
  int64_t nsites_sublat = nsites / n_sublat;
  int64_t n_leading = nsites_sublat;
  int64_t n_trailing = (n_sublat - 1) * nsites_sublat;

  for (auto prefix : combinatorics::Subsets<bit_t>(n_leading)) {

    int64_t nup_prefix = bits::popcount(prefix);
    int64_t nup_postfix = nup - nup_prefix;
    if ((nup_postfix < 0) || (nup_postfix > n_trailing)) {
      continue;
    }

    // if prefix is not rep, the full state also cannot be rep
    auto prefix_rep = action.reps_[n_sublat - 1][(int64_t)prefix];
    if (prefix_rep < prefix) {
      continue;
    }

#ifndef _OPENMP
    for (auto postfix :
         combinatorics::Combinations<bit_t>(n_trailing, nup_postfix)) {
      bit_t state = (prefix << n_trailing) | postfix;
      bit_t rep = action.representative(state);
      if (state == rep) {
        double norm = symmetries::norm(rep, action, characters);
        if (std::abs(norm) > 1e-6) {
          reps.push_back(rep);
          norms.push_back(norm);
        }
      }
    }
#else
    std::vector<std::vector<bit_t>> reps_thread;
    std::vector<std::vector<double>> norms_thread;
    auto combs = combinatorics::Combinations<bit_t>(n_trailing, nup_postfix);
    int64_t sz = combs.size();

#pragma omp parallel
    {
      int myid = omp_get_thread_num();
      int nthr = omp_get_num_threads();
#pragma omp single
      {
        reps_thread.resize(nthr);
        norms_thread.resize(nthr);
      }
      int64_t start = myid * (sz / nthr);
      int64_t end = (myid == nthr - 1) ? sz : (myid + 1) * (sz / nthr);
      for (auto it = combs.begin() + start, eit = combs.begin() + end;
           it != eit; ++it) {
        bit_t postfix = *it;
        bit_t state = (prefix << n_trailing) | postfix;
        bit_t rep = action.representative(state);
        if (state == rep) {
          double norm = symmetries::norm(rep, action, characters);
          if (std::abs(norm) > 1e-6) {
            reps_thread[myid].push_back(rep);
            norms_thread[myid].push_back(norm);
          }
        }
      }
    } // omp parallel

    for (auto const &rv : reps_thread)
      reps.insert(reps.end(), rv.begin(), rv.end());
    for (auto const &nv : norms_thread)
      norms.insert(norms.end(), nv.begin(), nv.end());
#endif
  } // for prefix

  reps.shrink_to_fit();
  norms.shrink_to_fit();
  return {reps, norms};
}

template <typename bit_t>
ska::flat_hash_map<bit_t, gsl::span<bit_t const>>
compute_rep_search_range(std::vector<bit_t> const &reps,
                         int64_t n_postfix_bits) {
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

} // namespace

template <typename bit_t, int n_sublat>
BasisSublattice<bit_t, n_sublat>::BasisSublattice(PermutationGroup const &group,
                                                  Vector const &characters) try
    : group_(group), action_(group_), characters_(characters),
      nsites_(group.nsites()), nup_(undefined),
      n_postfix_bits_(nsites_ - std::min(maximum_prefix_bits, nsites_)) {
  check_nsites_work_with_bits<bit_t>(nsites_);
  if (isreal(characters)) {
    std::tie(reps_, norms_) =
        reps_norms_no_sz(action_, characters.as<arma::vec>());
  } else {
    std::tie(reps_, norms_) =
        reps_norms_no_sz(action_, characters.as<arma::cx_vec>());
  }
}
XDIAG_CATCH

template <typename bit_t, int n_sublat>
BasisSublattice<bit_t, n_sublat>::BasisSublattice(int64_t nup,
                                                  PermutationGroup const &group,
                                                  Vector const &characters) try
    : group_(group), action_(group_), characters_(characters),
      nsites_(group.nsites()), nup_(nup),
      n_postfix_bits_(nsites_ - std::min(maximum_prefix_bits, nsites_)) {
  check_nsites_work_with_bits<bit_t>(nsites_);
  if (isreal(characters)) {
    std::tie(reps_, norms_) =
        reps_norms_sz(nup, action_, characters.as<arma::vec>());
  } else {
    std::tie(reps_, norms_) =
        reps_norms_sz(nup, action_, characters.as<arma::cx_vec>());
  }
}
XDIAG_CATCH

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
int64_t BasisSublattice<bit_t, n_sublat>::nsites() const {
  return nsites_;
}

template <typename bit_t, int n_sublat>
int64_t BasisSublattice<bit_t, n_sublat>::size() const {
  return (int64_t)reps_.size();
}

template <typename bit_t, int n_sublat>
int64_t BasisSublattice<bit_t, n_sublat>::d() const {
  return 2;
}

template <typename bit_t, int n_sublat>
Vector const &BasisSublattice<bit_t, n_sublat>::characters() const {
  return characters_;
}
template <typename bit_t, int n_sublat>

symmetries::SitePermutationSublattice<bit_t, n_sublat> const &
BasisSublattice<bit_t, n_sublat>::action() const {
  return action_;
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
  auto [rep, sym] = action_.representative_sym(state);
  int64_t idx = index_of_representative(rep);
  return {idx, sym};
}

template <typename bit_t, int n_sublat>
std::pair<int64_t, gsl::span<int64_t const>>
BasisSublattice<bit_t, n_sublat>::index_syms(bit_t state) const {
  auto [rep, syms] = action_.representative_syms(state);
  int64_t idx = index_of_representative(rep);
  return {idx, syms};
}

template <typename bit_t, int n_sublat>
bool BasisSublattice<bit_t, n_sublat>::operator==(
    BasisSublattice<bit_t, n_sublat> const &rhs) const {
  return (nsites_ == rhs.nsites_) && (nup_ == rhs.nup_) &&
         (n_postfix_bits_ == rhs.n_postfix_bits_) && (action_ == rhs.action_);
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

} // namespace xdiag::basis
