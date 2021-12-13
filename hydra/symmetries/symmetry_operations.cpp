#include "symmetry_operations.h"

#include <algorithm>
#include <fstream>

#include <hydra/common.h>

#include <hydra/utils/openmp_utils.h>

#include <hydra/combinatorics/bit_patterns.h>
#include <hydra/combinatorics/combinations_index.h>

#include <hydra/indexing/lintable.h>

#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/permutation_group_lookup.h>

#include <omp.h>
#include <unistd.h>

namespace hydra {

std::vector<std::vector<int>> read_permutations(std::string filename) {
  std::vector<std::vector<int>> lattice_symmetries;
  std::ifstream File(filename.c_str());

  if (File.fail()) {
    lila::Log.err("Error in read_spacegroup: Could not open file {}", filename);
    exit(EXIT_FAILURE);
  }

  std::string tobeparsed;
  std::string::size_type pos;

  // Jump to Sites and parse n_sites
  File >> tobeparsed;
  while (tobeparsed.find("[Sites]") == std::string::npos)
    File >> tobeparsed;
  pos = tobeparsed.find('=');
  int n_sites;
  if (pos != std::string::npos)
    n_sites = atoi(tobeparsed.substr(pos + 1, std::string::npos).c_str());
  else
    n_sites = -1;

  // Jump to SymmetryOps
  File >> tobeparsed;
  while (tobeparsed.find("[SymmetryOps]") == std::string::npos)
    File >> tobeparsed;

  // Read all symmetries
  int n_symmetries;
  pos = tobeparsed.find('=');
  if (pos != std::string::npos)
    n_symmetries = atoi(tobeparsed.substr(pos + 1, std::string::npos).c_str());
  else
    n_symmetries = -1;

  lattice_symmetries.resize(n_symmetries);
  for (int i = 0; i < n_symmetries; ++i) {
    File >> tobeparsed;
    for (int si = 0; si < n_sites; ++si) {
      int tosite;
      File >> tosite;
      lattice_symmetries[i].push_back(tosite);
    }
  }
  return lattice_symmetries;
}

} // namespace hydra

namespace hydra::symmetries {

bool is_valid_permutation(int n_sites, const int *permutation) {
  for (int i = 0; i < n_sites; ++i) {
    if (std::find(permutation, permutation + n_sites, i) ==
        permutation + n_sites)
      return false;
  }
  return true;
}

//////////////////////////////////////////////////////////
template <typename bit_t>
bit_t apply_permutation(bit_t state, int n_sites, const int *permutation) {
  bit_t tstate = 0;
  for (int site = 0; site < n_sites; ++site) {
    tstate |= ((state >> site) & 1) << permutation[site];
  }
  return tstate;
}

template uint16_t apply_permutation<uint16_t>(uint16_t, int, const int *);
template uint32_t apply_permutation<uint32_t>(uint32_t, int, const int *);
template uint64_t apply_permutation<uint64_t>(uint64_t, int, const int *);

//////////////////////////////////////////////////////////
template <typename bit_t, class GroupAction>
std::vector<int> stabilizer_symmetries(bit_t bits, GroupAction const &group) {
  std::vector<int> stable_syms;
  for (int sym = 0; sym < group.n_symmetries(); ++sym)
    if (group.apply(sym, bits) == bits)
      stable_syms.push_back(sym);
  return stable_syms;
}

template std::vector<int>
stabilizer_symmetries<uint16_t, PermutationGroupAction>(
    uint16_t bits, PermutationGroupAction const &group);
template std::vector<int>
stabilizer_symmetries<uint32_t, PermutationGroupAction>(
    uint32_t bits, PermutationGroupAction const &group);
template std::vector<int>
stabilizer_symmetries<uint64_t, PermutationGroupAction>(
    uint64_t bits, PermutationGroupAction const &group);
template std::vector<int>
stabilizer_symmetries<uint16_t, PermutationGroupLookup<uint16_t>>(
    uint16_t bits, PermutationGroupLookup<uint16_t> const &group);
template std::vector<int>
stabilizer_symmetries<uint32_t, PermutationGroupLookup<uint32_t>>(
    uint32_t bits, PermutationGroupLookup<uint32_t> const &group);
template std::vector<int>
stabilizer_symmetries<uint64_t, PermutationGroupLookup<uint64_t>>(
    uint64_t bits, PermutationGroupLookup<uint64_t> const &group);

//////////////////////////////////////////////////////////
template <typename bit_t, class GroupAction>
bit_t representative(bit_t state, GroupAction const &group_action) {
  if (group_action.n_symmetries() == 0)
    return state;

  bit_t rep = std::numeric_limits<bit_t>::max();
  for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < rep) {
      rep = tstate;
    }
  }
  return rep;
}

template uint16_t representative<uint16_t, PermutationGroupAction>(
    uint16_t state, PermutationGroupAction const &);
template uint32_t representative<uint32_t, PermutationGroupAction>(
    uint32_t state, PermutationGroupAction const &);
template uint64_t representative<uint64_t, PermutationGroupAction>(
    uint64_t state, PermutationGroupAction const &);
template uint16_t representative<uint16_t, PermutationGroupLookup<uint16_t>>(
    uint16_t state, PermutationGroupLookup<uint16_t> const &);
template uint32_t representative<uint32_t, PermutationGroupLookup<uint32_t>>(
    uint32_t state, PermutationGroupLookup<uint32_t> const &);
template uint64_t representative<uint64_t, PermutationGroupLookup<uint64_t>>(
    uint64_t state, PermutationGroupLookup<uint64_t> const &);

//////////////////////////////////////////////////////////
template <typename bit_t, class GroupAction>
bit_t representative_subset(bit_t state, GroupAction const &group_action,
                            gsl::span<int const> syms) {
  if (syms.size() == 0)
    return state;
  bit_t rep = std::numeric_limits<bit_t>::max();
  for (int sym : syms) {
    assert(sym < group_action.n_symmetries());
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < rep) {
      rep = tstate;
    }
  }
  return rep;
}

template uint16_t representative_subset<uint16_t, PermutationGroupAction>(
    uint16_t state, PermutationGroupAction const &, gsl::span<int const> syms);
template uint32_t representative_subset<uint32_t, PermutationGroupAction>(
    uint32_t state, PermutationGroupAction const &, gsl::span<int const> syms);
template uint64_t representative_subset<uint64_t, PermutationGroupAction>(
    uint64_t state, PermutationGroupAction const &, gsl::span<int const> syms);
template uint16_t
representative_subset<uint16_t, PermutationGroupLookup<uint16_t>>(
    uint16_t state, PermutationGroupLookup<uint16_t> const &,
    gsl::span<int const> syms);
template uint32_t
representative_subset<uint32_t, PermutationGroupLookup<uint32_t>>(
    uint32_t state, PermutationGroupLookup<uint32_t> const &,
    gsl::span<int const> syms);
template uint64_t
representative_subset<uint64_t, PermutationGroupLookup<uint64_t>>(
    uint64_t state, PermutationGroupLookup<uint64_t> const &,
    gsl::span<int const> syms);

//////////////////////////////////////////////////////////
template <typename bit_t, class GroupAction>
std::pair<bit_t, int> representative_sym(bit_t state,
                                         GroupAction const &group_action) {
  if (group_action.n_symmetries() == 0)
    return {state, 0};
  bit_t rep = std::numeric_limits<bit_t>::max();
  int rep_sym = 0;
  for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < rep) {
      rep = tstate;
      rep_sym = sym;
    }
  }
  return {rep, rep_sym};
}

template std::pair<uint16_t, int>
representative_sym<uint16_t, PermutationGroupAction>(
    uint16_t state, PermutationGroupAction const &);
template std::pair<uint32_t, int>
representative_sym<uint32_t, PermutationGroupAction>(
    uint32_t state, PermutationGroupAction const &);
template std::pair<uint64_t, int>
representative_sym<uint64_t, PermutationGroupAction>(
    uint64_t state, PermutationGroupAction const &);
template std::pair<uint16_t, int>
representative_sym<uint16_t, PermutationGroupLookup<uint16_t>>(
    uint16_t state, PermutationGroupLookup<uint16_t> const &);
template std::pair<uint32_t, int>
representative_sym<uint32_t, PermutationGroupLookup<uint32_t>>(
    uint32_t state, PermutationGroupLookup<uint32_t> const &);
template std::pair<uint64_t, int>
representative_sym<uint64_t, PermutationGroupLookup<uint64_t>>(
    uint64_t state, PermutationGroupLookup<uint64_t> const &);

//////////////////////////////////////////////////////////
template <typename bit_t, class GroupAction>
std::pair<bit_t, int> representative_sym_subset(bit_t state,
                                                GroupAction const &group_action,
                                                gsl::span<int const> syms) {
  if (syms.size() == 0)
    return {state, 0};
  bit_t rep = std::numeric_limits<bit_t>::max();
  int rep_sym = 0;
  for (auto sym : syms) {
    assert(sym < group_action.n_symmetries());
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < rep) {
      rep = tstate;
      rep_sym = sym;
    }
  }
  return {rep, rep_sym};
}

template std::pair<uint16_t, int>
representative_sym_subset<uint16_t, PermutationGroupAction>(
    uint16_t state, PermutationGroupAction const &, gsl::span<int const> syms);
template std::pair<uint32_t, int>
representative_sym_subset<uint32_t, PermutationGroupAction>(
    uint32_t state, PermutationGroupAction const &, gsl::span<int const> syms);
template std::pair<uint64_t, int>
representative_sym_subset<uint64_t, PermutationGroupAction>(
    uint64_t state, PermutationGroupAction const &, gsl::span<int const> syms);
template std::pair<uint16_t, int>
representative_sym_subset<uint16_t, PermutationGroupLookup<uint16_t>>(
    uint16_t state, PermutationGroupLookup<uint16_t> const &,
    gsl::span<int const> syms);
template std::pair<uint32_t, int>
representative_sym_subset<uint32_t, PermutationGroupLookup<uint32_t>>(
    uint32_t state, PermutationGroupLookup<uint32_t> const &,
    gsl::span<int const> syms);
template std::pair<uint64_t, int>
representative_sym_subset<uint64_t, PermutationGroupLookup<uint64_t>>(
    uint64_t state, PermutationGroupLookup<uint64_t> const &,
    gsl::span<int const> syms);

//////////////////////////////////////////////////////////
template <typename bit_t, class GroupAction>
std::pair<bit_t, std::vector<int>>
representative_syms(bit_t state, GroupAction const &group_action) {
  if (group_action.n_symmetries() == 0)
    return {state, std::vector<int>()};

  bit_t rep = representative(state, group_action);
  std::vector<int> rep_syms;
  for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate == rep) {
      rep_syms.push_back(sym);
    }
  }
  return {rep, rep_syms};
}

template std::pair<uint16_t, std::vector<int>>
representative_syms<uint16_t, PermutationGroupAction>(
    uint16_t state, PermutationGroupAction const &);
template std::pair<uint32_t, std::vector<int>>
representative_syms<uint32_t, PermutationGroupAction>(
    uint32_t state, PermutationGroupAction const &);
template std::pair<uint64_t, std::vector<int>>
representative_syms<uint64_t, PermutationGroupAction>(
    uint64_t state, PermutationGroupAction const &);
template std::pair<uint16_t, std::vector<int>>
representative_syms<uint16_t, PermutationGroupLookup<uint16_t>>(
    uint16_t state, PermutationGroupLookup<uint16_t> const &);
template std::pair<uint32_t, std::vector<int>>
representative_syms<uint32_t, PermutationGroupLookup<uint32_t>>(
    uint32_t state, PermutationGroupLookup<uint32_t> const &);
template std::pair<uint64_t, std::vector<int>>
representative_syms<uint64_t, PermutationGroupLookup<uint64_t>>(
    uint64_t state, PermutationGroupLookup<uint64_t> const &);

//////////////////////////////////////////////////////////
template <typename bit_t, class GroupAction>
std::vector<int> mapping_syms(bit_t origin, bit_t target,
                              GroupAction const &group_action) {
  std::vector<int> syms;
  for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, origin);
    if (tstate == target) {
      syms.push_back(sym);
    }
  }
  return syms;
}

template std::vector<int>
mapping_syms<uint16_t, PermutationGroupAction>(uint16_t, uint16_t,
                                               PermutationGroupAction const &);
template std::vector<int>
mapping_syms<uint32_t, PermutationGroupAction>(uint32_t, uint32_t,
                                               PermutationGroupAction const &);
template std::vector<int>
mapping_syms<uint64_t, PermutationGroupAction>(uint64_t, uint64_t,
                                               PermutationGroupAction const &);
template std::vector<int>
mapping_syms<uint16_t, PermutationGroupLookup<uint16_t>>(
    uint16_t, uint16_t, PermutationGroupLookup<uint16_t> const &);
template std::vector<int>
mapping_syms<uint32_t, PermutationGroupLookup<uint32_t>>(
    uint32_t, uint32_t, PermutationGroupLookup<uint32_t> const &);
template std::vector<int>
mapping_syms<uint64_t, PermutationGroupLookup<uint64_t>>(
    uint64_t, uint64_t, PermutationGroupLookup<uint64_t> const &);

//////////////////////////////////////////////////////////
template <typename bit_t, class GroupAction>
double norm(bit_t state, GroupAction const &group_action,
            Representation const &irrep) {
  complex amplitude = 0.0;
  for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate == state) {
      amplitude += irrep.character(sym);
    }
  }
  return std::sqrt(std::abs(amplitude));
}

template double
norm<uint16_t, PermutationGroupAction>(uint16_t, PermutationGroupAction const &,
                                       Representation const &);
template double
norm<uint32_t, PermutationGroupAction>(uint32_t, PermutationGroupAction const &,
                                       Representation const &);
template double
norm<uint64_t, PermutationGroupAction>(uint64_t, PermutationGroupAction const &,
                                       Representation const &);
template double norm<uint16_t, PermutationGroupLookup<uint16_t>>(
    uint16_t, PermutationGroupLookup<uint16_t> const &, Representation const &);
template double norm<uint32_t, PermutationGroupLookup<uint32_t>>(
    uint32_t, PermutationGroupLookup<uint32_t> const &, Representation const &);
template double norm<uint64_t, PermutationGroupLookup<uint64_t>>(
    uint64_t, PermutationGroupLookup<uint64_t> const &, Representation const &);

//////////////////////////////////////////////////////////
template <typename bit_t, class GroupAction>
double norm_fermionic(bit_t state, GroupAction const &group_action,
                      Representation const &irrep) {
  complex amplitude = 0.0;
  int n_sites = group_action.n_sites();
  auto work = fermi_work(n_sites);
  const int *sym_ptr = group_action.permutation_array().data();
  for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate == state) {
      if (fermi_bool_of_permutation(state, sym_ptr, work.data())) {
        amplitude -= irrep.character(sym);
      } else {
        amplitude += irrep.character(sym);
      }
    }
    sym_ptr += n_sites;
  }
  return std::sqrt(std::abs(amplitude));
}

template double norm_fermionic<uint16_t, PermutationGroupAction>(
    uint16_t state, PermutationGroupAction const &, Representation const &);
template double norm_fermionic<uint32_t, PermutationGroupAction>(
    uint32_t state, PermutationGroupAction const &, Representation const &);
template double norm_fermionic<uint64_t, PermutationGroupAction>(
    uint64_t state, PermutationGroupAction const &, Representation const &);
template double norm_fermionic<uint16_t, PermutationGroupLookup<uint16_t>>(
    uint16_t state, PermutationGroupLookup<uint16_t> const &,
    Representation const &);
template double norm_fermionic<uint32_t, PermutationGroupLookup<uint32_t>>(
    uint32_t state, PermutationGroupLookup<uint32_t> const &,
    Representation const &);
template double norm_fermionic<uint64_t, PermutationGroupLookup<uint64_t>>(
    uint64_t state, PermutationGroupLookup<uint64_t> const &,
    Representation const &);

//////////////////////////////////////////////////////////
template <typename bit_t, class GroupAction>
double norm_electron(bit_t ups, bit_t dns, GroupAction const &group_action,
                     Representation const &irrep) {
  assert(group_action.n_symmetries() == irrep.size());
  complex amplitude = 0.0;
  int n_sites = group_action.n_sites();
  auto work = fermi_work(n_sites);
  const int *sym_ptr = group_action.permutation_array().data();
  for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tups = group_action.apply(sym, ups);

    if (tups == ups) {
      bool fermi_bool_ups =
          fermi_bool_of_permutation(ups, sym_ptr, work.data());

      bit_t tdns = group_action.apply(sym, dns);

      if (tdns == dns) {
        bool fermi_bool_dns =
            fermi_bool_of_permutation(dns, sym_ptr, work.data());

        if (fermi_bool_ups == fermi_bool_dns) {
          amplitude += irrep.character(sym);
        } else {
          amplitude -= irrep.character(sym);
        }
      }
    }
    sym_ptr += n_sites;
  }
  return std::sqrt(std::abs(amplitude));
}

template double norm_electron<uint16_t, PermutationGroupAction>(
    uint16_t, uint16_t, PermutationGroupAction const &, Representation const &);
template double norm_electron<uint32_t, PermutationGroupAction>(
    uint32_t, uint32_t, PermutationGroupAction const &, Representation const &);
template double norm_electron<uint64_t, PermutationGroupAction>(
    uint64_t, uint64_t, PermutationGroupAction const &, Representation const &);
template double norm_electron<uint16_t, PermutationGroupLookup<uint16_t>>(
    uint16_t, uint16_t, PermutationGroupLookup<uint16_t> const &,
    Representation const &);
template double norm_electron<uint32_t, PermutationGroupLookup<uint32_t>>(
    uint32_t, uint32_t, PermutationGroupLookup<uint32_t> const &,
    Representation const &);
template double norm_electron<uint64_t, PermutationGroupLookup<uint64_t>>(
    uint64_t, uint64_t, PermutationGroupLookup<uint64_t> const &,
    Representation const &);

//////////////////////////////////////////////////////////
template <typename bit_t, class GroupAction>
double
norm_electron_subset(bit_t ups, bit_t dns, GroupAction const &group_action,
                     Representation const &irrep, gsl::span<int const> syms) {
  assert(group_action.n_symmetries() == irrep.size());
  complex amplitude = 0.0;
  int n_sites = group_action.n_sites();
  auto work = fermi_work(n_sites);
  const int *sym_ptr = group_action.permutation_array().data();
  for (int sym : syms) {
    assert(sym < group_action.n_symmetries());

    const int *sym_ptr_current = sym_ptr + sym * n_sites;
    bit_t tups = group_action.apply(sym, ups);

    if (tups == ups) {
      bool fermi_bool_ups =
          fermi_bool_of_permutation(ups, sym_ptr_current, work.data());

      bit_t tdns = group_action.apply(sym, dns);

      if (tdns == dns) {
        bool fermi_bool_dns =
            fermi_bool_of_permutation(dns, sym_ptr_current, work.data());
        if (fermi_bool_ups == fermi_bool_dns) {
          amplitude += irrep.character(sym);
        } else {
          amplitude -= irrep.character(sym);
        }
      }
    }
  }
  return std::sqrt(std::abs(amplitude));
}

template double norm_electron_subset<uint16_t, PermutationGroupAction>(
    uint16_t, uint16_t, PermutationGroupAction const &, Representation const &,
    gsl::span<int const>);
template double norm_electron_subset<uint32_t, PermutationGroupAction>(
    uint32_t, uint32_t, PermutationGroupAction const &, Representation const &,
    gsl::span<int const>);
template double norm_electron_subset<uint64_t, PermutationGroupAction>(
    uint64_t, uint64_t, PermutationGroupAction const &, Representation const &,
    gsl::span<int const>);
template double
norm_electron_subset<uint16_t, PermutationGroupLookup<uint16_t>>(
    uint16_t, uint16_t, PermutationGroupLookup<uint16_t> const &,
    Representation const &, gsl::span<int const>);
template double
norm_electron_subset<uint32_t, PermutationGroupLookup<uint32_t>>(
    uint32_t, uint32_t, PermutationGroupLookup<uint32_t> const &,
    Representation const &, gsl::span<int const>);
template double
norm_electron_subset<uint64_t, PermutationGroupLookup<uint64_t>>(
    uint64_t, uint64_t, PermutationGroupLookup<uint64_t> const &,
    Representation const &, gsl::span<int const>);

//////////////////////////////////////////////////////////
template <typename bit_t, class GroupAction, class LinTable>
std::tuple<std::vector<bit_t>, std::vector<idx_t>, std::vector<int>,
           std::vector<std::pair<span_size_t, span_size_t>>>
representatives_indices_symmetries_limits(int n_par,
                                          GroupAction const &group_action,
                                          LinTable const &lintable) {

  using combinatorics::binomial;
  using combinatorics::CombinationsIndex;
  using combinatorics::CombinationsIndexThread;

  int n_sites = group_action.n_sites();
  std::vector<idx_t> idces(binomial(n_sites, n_par), invalid_index);
  std::vector<std::pair<span_size_t, span_size_t>> sym_limits(
      binomial(n_sites, n_par));

  assert(n_par <= n_sites);
  assert(n_par >= 0);

#ifndef _OPENMP

  // Compute all representatives
  std::vector<bit_t> reps;
  for (auto [state, idx] : CombinationsIndex<bit_t>(n_sites, n_par)) {
    bit_t rep = representative(state, group_action);
    if (rep == state) {
      idces[idx] = reps.size();
      reps.push_back(rep);
    }
  }

  // Compute indices of up-representatives and stabilizer symmetries
  std::vector<int> syms;
  for (auto [state, idx] : CombinationsIndex<bit_t>(n_sites, n_par)) {
    bit_t rep = representative(state, group_action);
    idces[idx] = idces[lintable.index(rep)];

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
#else

  // Compute all representatives
  std::vector<std::vector<bit_t>> reps_thread;
  std::vector<std::vector<int>> syms_thread;

#pragma omp parallel
  {
    int myid = omp_get_thread_num();
    int rank = omp_get_num_threads();

#pragma omp single
    { reps_thread.resize(rank); }

    for (auto [state, idx] : CombinationsIndexThread<bit_t>(n_sites, n_par)) {
      bit_t rep = representative(state, group_action);
      if (rep == state) {
        idces[idx] = reps_thread[myid].size(); // just local index
        reps_thread[myid].push_back(rep);
      }
    }

#pragma omp barrier

    idx_t offset = 0;
    for (int id = 0; id < myid; ++id) {
      offset += reps_thread[id].size();
    }

    for (auto [state, idx] : CombinationsIndexThread<bit_t>(n_sites, n_par)) {
      if (idces[idx] != invalid_index) {
        idces[idx] += offset;
      }
    }
  } // pragma omp parallel

  std::vector<bit_t> reps = utils::combine_vectors(reps_thread);

  // Compute indices of up-representatives and stabilizer symmetries
#pragma omp parallel
  {
    int myid = omp_get_thread_num();
    int rank = omp_get_num_threads();

#pragma omp single
    { syms_thread.resize(rank); }

    for (auto [state, idx] : CombinationsIndexThread<bit_t>(n_sites, n_par)) {
      bit_t rep = representative(state, group_action);

      assert(bitops::popcnt(rep) == n_par);
      assert(rep < (bit_t)1 << n_sites);
      assert(rep <= state);

      idces[idx] = idces[lintable.index(rep)];

      assert(idces[idx] != invalid_index);

      // Determine the symmetries that yield the up-representative
      span_size_t begin = syms_thread[myid].size();
      for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
        if (group_action.apply(sym, state) == rep)
          syms_thread[myid].push_back(sym);
      }
      span_size_t end = syms_thread[myid].size();
      sym_limits[idx] = {begin, end - begin};
    }

#pragma omp barrier

    span_size_t offset = 0;
    for (int id = 0; id < myid; ++id) {
      offset += syms_thread[id].size();
    }

    for (auto [state, idx] : CombinationsIndexThread<bit_t>(n_sites, n_par)) {
      auto [begin, length] = sym_limits[idx];
      sym_limits[idx] = {begin + offset, length};
    }

  } // pragma omp parallel

  // Combine syms
  std::vector<int> syms = utils::combine_vectors(syms_thread);

#endif

  // lila::Log("reps.size(): {}, syms.size(): {}", reps.size(), syms.size());
  // for (auto rep : reps)
  //   lila::Log("rep: {}", rep);
  // for (auto idx : idces)
  //   lila::Log("idx: {}", idx);

  // for (auto sym : syms)
  //   lila::Log("sym: {}", sym);

  // for (auto [b, l] : sym_limits) {

  //   lila::Log("lim: {} {}", b, l);
  // }
  return {reps, idces, syms, sym_limits};
}

template std::tuple<std::vector<uint16_t>, std::vector<idx_t>, std::vector<int>,
                    std::vector<std::pair<span_size_t, span_size_t>>>
representatives_indices_symmetries_limits<uint16_t, PermutationGroupAction,
                                          indexing::LinTable<uint16_t>>(
    int, PermutationGroupAction const &, indexing::LinTable<uint16_t> const &);

template std::tuple<std::vector<uint32_t>, std::vector<idx_t>, std::vector<int>,
                    std::vector<std::pair<span_size_t, span_size_t>>>
representatives_indices_symmetries_limits<uint32_t, PermutationGroupAction,
                                          indexing::LinTable<uint32_t>>(
    int, PermutationGroupAction const &, indexing::LinTable<uint32_t> const &);

template std::tuple<std::vector<uint64_t>, std::vector<idx_t>, std::vector<int>,
                    std::vector<std::pair<span_size_t, span_size_t>>>
representatives_indices_symmetries_limits<uint64_t, PermutationGroupAction,
                                          indexing::LinTable<uint64_t>>(
    int, PermutationGroupAction const &, indexing::LinTable<uint64_t> const &);
template std::tuple<std::vector<uint16_t>, std::vector<idx_t>, std::vector<int>,
                    std::vector<std::pair<span_size_t, span_size_t>>>
representatives_indices_symmetries_limits<
    uint16_t, PermutationGroupLookup<uint16_t>, indexing::LinTable<uint16_t>>(
    int, PermutationGroupLookup<uint16_t> const &,
    indexing::LinTable<uint16_t> const &);

template std::tuple<std::vector<uint32_t>, std::vector<idx_t>, std::vector<int>,
                    std::vector<std::pair<span_size_t, span_size_t>>>
representatives_indices_symmetries_limits<
    uint32_t, PermutationGroupLookup<uint32_t>, indexing::LinTable<uint32_t>>(
    int, PermutationGroupLookup<uint32_t> const &,
    indexing::LinTable<uint32_t> const &);

template std::tuple<std::vector<uint64_t>, std::vector<idx_t>, std::vector<int>,
                    std::vector<std::pair<span_size_t, span_size_t>>>
representatives_indices_symmetries_limits<
    uint64_t, PermutationGroupLookup<uint64_t>, indexing::LinTable<uint64_t>>(
    int, PermutationGroupLookup<uint64_t> const &,
    indexing::LinTable<uint64_t> const &);

//////////////////////////////////////////////////////////
template <typename bit_t, class GroupAction>
std::vector<bool> fermi_bool_table(int npar, GroupAction const &group_action) {

  using combinatorics::binomial;
  using combinatorics::CombinationsIndex;
  using combinatorics::CombinationsIndexThread;
  using combinatorics::get_next_pattern;
  using combinatorics::get_nth_pattern;

  int n_sites = group_action.n_sites();
  int n_symmetries = group_action.n_symmetries();
  idx_t raw_size = binomial(n_sites, npar);
  std::vector<bool> fermi_bool_table(raw_size * n_symmetries);
  const int *sym_ptr = group_action.permutation_array().data();

#ifndef _OPENMP
  auto fermi_work = symmetries::fermi_work(n_sites);
  for (int sym = 0; sym < n_symmetries; ++sym) {
    for (auto [state, idx] : CombinationsIndex<bit_t>(n_sites, npar)) {
      fermi_bool_table[sym * raw_size + idx] =
          fermi_bool_of_permutation(state, sym_ptr, fermi_work.data());
    }
    sym_ptr += n_sites;
  }
#else

  for (int sym = 0; sym < n_symmetries; ++sym) {

    std::vector<std::vector<bool>> fermi_bool_table_local;

    // auto t1 = lila::rightnow();
#pragma omp parallel
    {
      int myid = omp_get_thread_num();
      int rank = omp_get_num_threads();

      auto fermi_work = symmetries::fermi_work(n_sites);

      auto comb = CombinationsIndexThread<bit_t>(n_sites, npar);
      auto begin = comb.begin();
      auto end = comb.end();
      idx_t length = end.idx() - begin.idx();

#pragma omp single
      { fermi_bool_table_local.resize(rank); }

#pragma omp barrier

      fermi_bool_table_local[myid].resize(length);
      auto [state, idxl] = *begin;

      for (idx_t idx = 0; idx < length; ++idx) {
        bool fermi_bool =
            fermi_bool_of_permutation(state, sym_ptr, fermi_work.data());
        fermi_bool_table_local[myid][idx] = fermi_bool;
        state = get_next_pattern(state);
      }
    } // pragma omp parallel
    // lila::timing(t1, lila::rightnow(), "fill");

    // auto t2 = lila::rightnow();
    auto fermi_bool_table_for_sym =
        utils::combine_vectors_copy(fermi_bool_table_local);
    std::copy(fermi_bool_table_for_sym.begin(), fermi_bool_table_for_sym.end(),
              fermi_bool_table.begin() + sym * raw_size);
    // lila::timing(t2, lila::rightnow(), "combine");
    sym_ptr += n_sites;
  }

#endif // ifndef _OPENMP

  return fermi_bool_table;
}

template std::vector<bool> fermi_bool_table<uint16_t, PermutationGroupAction>(
    int, PermutationGroupAction const &);
template std::vector<bool> fermi_bool_table<uint32_t, PermutationGroupAction>(
    int, PermutationGroupAction const &);
template std::vector<bool> fermi_bool_table<uint64_t, PermutationGroupAction>(
    int, PermutationGroupAction const &);
template std::vector<bool>
fermi_bool_table<uint16_t, PermutationGroupLookup<uint16_t>>(
    int, PermutationGroupLookup<uint16_t> const &);
template std::vector<bool>
fermi_bool_table<uint32_t, PermutationGroupLookup<uint32_t>>(
    int, PermutationGroupLookup<uint32_t> const &);
template std::vector<bool>
fermi_bool_table<uint64_t, PermutationGroupLookup<uint64_t>>(
    int, PermutationGroupLookup<uint64_t> const &);

} // namespace hydra::symmetries
