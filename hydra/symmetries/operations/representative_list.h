#pragma once

#include <utility>
#include <vector>

#include <lila/external/gsl/span>

#include <hydra/common.h>
#include <hydra/symmetries/operations/group_action_operations.h>
#include <hydra/symmetries/operations/symmetry_operations.h>
#include <hydra/utils/logger.h>

namespace hydra::symmetries {

using span_size_t = gsl::span<int const>::size_type;

template <typename bit_t, class StatesIndexing, class GroupAction>
inline std::tuple<std::vector<bit_t>, std::vector<idx_t>, std::vector<int>,
                  std::vector<std::pair<span_size_t, span_size_t>>>
representatives_indices_symmetries_limits(StatesIndexing &&states_indexing,
                                          GroupAction &&group_action) {

  using combinatorics::binomial;

  idx_t size = states_indexing.size();
  std::vector<idx_t> idces(size);
  std::vector<std::pair<span_size_t, span_size_t>> sym_limits(size);

  // Compute all representatives
  std::vector<bit_t> reps;
  for (auto [state, idx] : states_indexing.states_indices()) {
    if (is_representative(state, group_action)) {
      idces[idx] = reps.size();
      reps.push_back(state);
    }
  }
  reps.shrink_to_fit();
  
  // Compute indices of up-representatives and stabilizer symmetries
  std::vector<int> syms;
  for (auto [state, idx] : states_indexing.states_indices()) {
    bit_t rep = representative(state, group_action);
    idces[idx] = idces[states_indexing.index(rep)];

    assert(idces[idx] != invalid_index);

    // Determine the symmetries that yield the up-representative
    span_size_t begin = syms.size();
    for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
      if (group_action.apply(sym, state) == rep)
        syms.push_back(sym);
    }
    span_size_t end = syms.size();
    sym_limits[idx] = {begin, end - begin};
  }

  //   // Compute indices of up-representatives and stabilizer symmetries
  // std::vector<int> syms;
  // for (idx_t rep_idx=0; rep_idx<reps.size(); ++rep_idx) {
  //   bit_t rep = reps[rep_idx];

  //   std::map<bit_t, std::vector<int>> orbits;
    
  //   for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
  //     bit_t state = group_action.apply(sym, rep);
  //     if (std::find(orbits.begin(), orbits.end(), state) == orbits.end()){
  // 	orbits[state] = {sym}
  //     } else {
  // 	orbits[state].push_back(sym);
  //     }

  //     idx_t idx = states_indexing.index(state);
  //     idces[idx] = idces[rep_idx];

      
  //   }
      
  //   bit_t rep = representative(state, group_action);
  //   idces[idx] = idces[states_indexing.index(rep)];

  //   assert(idces[idx] != invalid_index);

  //   // Determine the symmetries that yield the up-representative
  //   span_size_t begin = syms.size();
  //   for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
  //     if (group_action.apply(sym, state) == rep)
  //       syms.push_back(sym);
  //   }
  //   span_size_t end = syms.size();
  //   sym_limits[idx] = {begin, end - begin};
  // }


  return {reps, idces, syms, sym_limits};
}

template <typename bit_t, class StatesIndexing, class GroupAction>
inline std::tuple<std::vector<bit_t>, std::vector<idx_t>, std::vector<int>,
                  std::vector<std::pair<span_size_t, span_size_t>>,
                  std::vector<double>>
representatives_indices_symmetries_limits_norms(
    StatesIndexing &&states_indexing, GroupAction &&group_action,
    Representation const &irrep) {
  idx_t size = states_indexing.size();

  std::vector<bit_t> reps;
  std::vector<idx_t> idces(size, invalid_index);
  std::vector<int> syms;
  std::vector<std::pair<span_size_t, span_size_t>> sym_limits(size);
  std::vector<double> norms;

  // Compute all representatives
  for (auto [state, idx] : states_indexing.states_indices()) {
    // register state if it's a representative and has non-zero norm
    if (is_representative(state, group_action)) {
      double norm = symmetries::norm(state, group_action, irrep);
      if (norm > 1e-6) {
        idces[idx] = reps.size();
        reps.push_back(state);
        norms.push_back(norm);
      }
    }
  }

  // Compute indices of up-representatives and stabilizer symmetries
  for (auto [state, idx] : states_indexing.states_indices()) {

    bit_t rep = representative(state, group_action);
    idx_t rep_idx = idces[states_indexing.index(rep)];
    if (rep_idx != invalid_index) { // can be invalid if zero norm
      idces[idx] = rep_idx;

      // Add syms yielding the representative
      std::vector<int> rep_syms = mapping_syms(state, rep, group_action);
      span_size_t start = syms.size();
      syms.insert(syms.end(), rep_syms.begin(), rep_syms.end());
      span_size_t end = syms.size();
      sym_limits[idx] = {start, end - start};
    }
  }

  // Resize to correct size
  reps.shrink_to_fit();
  idces.shrink_to_fit();
  syms.shrink_to_fit();
  sym_limits.shrink_to_fit();
  norms.shrink_to_fit();

  return {reps, idces, syms, sym_limits, norms};
}

template <typename bit_t, class States, class GroupAction>
inline std::tuple<std::vector<bit_t>, std::vector<double>,
                  std::vector<std::pair<span_size_t, span_size_t>>,
                  std::vector<idx_t>, idx_t>
electron_dns_norms_limits_offset_size(std::vector<bit_t> const &reps_up,
                                      States &&states_dns,
                                      GroupAction &&group_action,
                                      Representation const &irrep) {

  std::vector<bit_t> dns_storage;
  std::vector<double> norms_storage;
  std::vector<std::pair<span_size_t, span_size_t>> dns_limits(reps_up.size());
  std::vector<idx_t> ups_offset((reps_up.size()));

  // if ups have trivial stabilizer, dns  are stored in front
  for (bit_t dns : states_dns) {
    dns_storage.push_back(dns);
    norms_storage.push_back(1.0);
  }

  idx_t size = 0;
  idx_t idx_up = 0;
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
} // namespace hydra::symmetries
