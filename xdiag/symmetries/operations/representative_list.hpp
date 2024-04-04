#pragma once

#include <numeric>
#include <utility>
#include <vector>

#include <xdiag/extern/gsl/span>
#include <xdiag/common.hpp>
#include <xdiag/symmetries/operations/group_action_operations.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::symmetries {

using span_size_t = gsl::span<int64_t const>::size_type;

template <typename bit_t, class StatesIndexing, class GroupAction>
inline std::tuple<
    std::vector<bit_t>, std::vector<int64_t>, std::vector<int64_t>,
    std::vector<std::pair<span_size_t, span_size_t>>, std::vector<double>>
representatives_indices_symmetries_limits_norms(
    StatesIndexing &&states_indexing, GroupAction &&group_action,
    Representation const &irrep) {
  int64_t size = states_indexing.size();
  std::vector<int64_t> idces(size, invalid_index);

  // Compute all representatives
  std::vector<bit_t> reps;
  std::vector<double> norms;
  for (auto [state, idx] : states_indexing.states_indices()) {
    if (is_representative(state, group_action)) {
      double nrm = symmetries::norm(state, group_action, irrep);
      if (std::abs(nrm) > 1e-6) {
        idces[idx] = reps.size();
        reps.push_back(state);
        norms.push_back(nrm);
      }
    }
  }
  reps.shrink_to_fit();
  norms.shrink_to_fit();
  int64_t n_reps = reps.size();

  // Determine the number of syms yielding the representative for each state
  std::vector<int64_t> n_syms_for_state(size, 0);
  for (int64_t rep_idx = 0; rep_idx < n_reps; ++rep_idx) {
    bit_t rep = reps[rep_idx];

    for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
      bit_t state = group_action.apply(sym, rep);
      int64_t idx = states_indexing.index(state);
      idces[idx] = rep_idx;
      ++n_syms_for_state[idx];
    }
  }

  // compute size and allocate syms array
  int64_t n_syms =
      std::accumulate(n_syms_for_state.begin(), n_syms_for_state.end(), 0);
  std::vector<int64_t> syms(n_syms, 0);

  // compute the sym offsets
  std::vector<int64_t> n_syms_for_state_offset(size, 0);

  // std::exclusive_scan(n_syms_for_state.begin(), n_syms_for_state.end(),
  //                     n_syms_for_state_offset.begin(), 0);

  std::partial_sum(n_syms_for_state.begin(), n_syms_for_state.end() - 1,
                   n_syms_for_state_offset.begin() + 1);

  // set the sym_limits
  std::vector<std::pair<span_size_t, span_size_t>> sym_limits(size);
  for (int64_t idx = 0; idx < size; ++idx) {
    sym_limits[idx] = {n_syms_for_state_offset[idx], n_syms_for_state[idx]};
  }

  // empty n_syms_for_state again, now used as a counter
  std::fill(n_syms_for_state.begin(), n_syms_for_state.end(), 0);

  // calculate the symmetries yielding the representative
  auto const &group = group_action.permutation_group();
  for (int64_t rep_idx = 0; rep_idx < n_reps; ++rep_idx) {
    bit_t rep = reps[rep_idx];

    for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
      bit_t state = group_action.apply(sym, rep);
      int64_t idx = states_indexing.index(state);

      int64_t sym_inv = group.inverse(sym);
      int64_t idx_sym = n_syms_for_state_offset[idx] + n_syms_for_state[idx]++;
      syms[idx_sym] = sym_inv;
    }
  }

  return {reps, idces, syms, sym_limits, norms};
}

template <typename bit_t, class StatesIndexing, class GroupAction>
inline std::tuple<std::vector<bit_t>, std::vector<int64_t>,
                  std::vector<int64_t>,
                  std::vector<std::pair<span_size_t, span_size_t>>>
representatives_indices_symmetries_limits(StatesIndexing &&states_indexing,
                                          GroupAction &&group_action) {
  auto irrep = trivial_representation(group_action.n_symmetries());
  auto [reps, idces, syms, sym_limits, norms] =
      representatives_indices_symmetries_limits_norms<bit_t>(
          states_indexing, group_action, irrep);
  (void)norms;
  return {reps, idces, syms, sym_limits};
}

template <typename bit_t, class States, class GroupAction>
inline std::tuple<std::vector<bit_t>, std::vector<double>,
                  std::vector<std::pair<span_size_t, span_size_t>>,
                  std::vector<int64_t>, int64_t>
electron_dns_norms_limits_offset_size(std::vector<bit_t> const &reps_up,
                                      States &&states_dns,
                                      GroupAction &&group_action,
                                      Representation const &irrep) {

  std::vector<bit_t> dns_storage;
  std::vector<double> norms_storage;
  std::vector<std::pair<span_size_t, span_size_t>> dns_limits(reps_up.size());
  std::vector<int64_t> ups_offset((reps_up.size()));

  // if ups have trivial stabilizer, dns  are stored in front
  for (bit_t dns : states_dns) {
    dns_storage.push_back(dns);
    norms_storage.push_back(1.0);
  }

  int64_t size = 0;
  int64_t idx_up = 0;
  for (bit_t ups : reps_up) {
    ups_offset[idx_up] = size;
    auto syms = mapping_syms(ups, ups, group_action);

    // ups have trivial stabilizer -> dns stored in beginning
    if (syms.size() == 1) {
      span_size_t start = 0;
      span_size_t length = states_dns.size();
      dns_limits[idx_up] = {start, length};
      size += length;
    } // ups have non-trivial stabilizer, we store the dns configurations
    else {
      span_size_t start = dns_storage.size();
      for (bit_t dns : states_dns) {
        bit_t dns_rep = representative_subset(dns, group_action, syms);
        if (dns == dns_rep) {
          double norm =
              norm_electron_subset(ups, dns, group_action, irrep, syms);
          if (norm > 1e-6) {
            dns_storage.push_back(dns_rep);
            norms_storage.push_back(norm);
          }
        }
      }
      span_size_t end = dns_storage.size();
      span_size_t length = end - start;
      dns_limits[idx_up] = {start, length};
      size += length;
    }
    ++idx_up;
  }

  return {dns_storage, norms_storage, dns_limits, ups_offset, size};
}
} // namespace xdiag::symmetries
