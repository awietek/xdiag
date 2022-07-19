#include "symmetry_operations.h"

#include <algorithm>
#include <fstream>

#include <hydra/common.h>

#include <hydra/utils/logger.h>
#include <hydra/utils/openmp_utils.h>

#include <hydra/combinatorics/bit_patterns.h>
#include <hydra/combinatorics/combinations_index.h>

#include <hydra/indexing/lintable.h>

#include <hydra/symmetries/group_action/group_action.h>
#include <hydra/symmetries/group_action/group_action_lookup.h>

namespace hydra {

std::vector<std::vector<int>> read_permutations(std::string filename) {
  std::vector<std::vector<int>> lattice_symmetries;
  std::ifstream File(filename.c_str());

  if (File.fail()) {
    Log.err("Error in read_spacegroup: Could not open file {}", filename);
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
template <typename bit_t>
bit_t apply_permutation(bit_t state, gsl::span<int const> permutation) {
  bit_t tstate = 0;
  for (int site = 0; site < permutation.size(); ++site) {
    tstate |= ((state >> site) & 1) << permutation[site];
  }
  return tstate;
}

template uint16_t apply_permutation<uint16_t>(uint16_t, gsl::span<int const>);
template uint32_t apply_permutation<uint32_t>(uint32_t, gsl::span<int const>);
template uint64_t apply_permutation<uint64_t>(uint64_t, gsl::span<int const>);

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
stabilizer_symmetries<uint16_t, GroupAction>(uint16_t bits,
                                             GroupAction const &group);
template std::vector<int>
stabilizer_symmetries<uint32_t, GroupAction>(uint32_t bits,
                                             GroupAction const &group);
template std::vector<int>
stabilizer_symmetries<uint64_t, GroupAction>(uint64_t bits,
                                             GroupAction const &group);
template std::vector<int>
stabilizer_symmetries<uint16_t, GroupActionLookup<uint16_t>>(
    uint16_t bits, GroupActionLookup<uint16_t> const &group);
template std::vector<int>
stabilizer_symmetries<uint32_t, GroupActionLookup<uint32_t>>(
    uint32_t bits, GroupActionLookup<uint32_t> const &group);
template std::vector<int>
stabilizer_symmetries<uint64_t, GroupActionLookup<uint64_t>>(
    uint64_t bits, GroupActionLookup<uint64_t> const &group);

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

template uint16_t representative<uint16_t, GroupAction>(uint16_t state,
                                                        GroupAction const &);
template uint32_t representative<uint32_t, GroupAction>(uint32_t state,
                                                        GroupAction const &);
template uint64_t representative<uint64_t, GroupAction>(uint64_t state,
                                                        GroupAction const &);
template uint16_t representative<uint16_t, GroupActionLookup<uint16_t>>(
    uint16_t state, GroupActionLookup<uint16_t> const &);
template uint32_t representative<uint32_t, GroupActionLookup<uint32_t>>(
    uint32_t state, GroupActionLookup<uint32_t> const &);
template uint64_t representative<uint64_t, GroupActionLookup<uint64_t>>(
    uint64_t state, GroupActionLookup<uint64_t> const &);

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

template uint16_t representative_subset<uint16_t, GroupAction>(
    uint16_t state, GroupAction const &, gsl::span<int const> syms);
template uint32_t representative_subset<uint32_t, GroupAction>(
    uint32_t state, GroupAction const &, gsl::span<int const> syms);
template uint64_t representative_subset<uint64_t, GroupAction>(
    uint64_t state, GroupAction const &, gsl::span<int const> syms);
template uint16_t representative_subset<uint16_t, GroupActionLookup<uint16_t>>(
    uint16_t state, GroupActionLookup<uint16_t> const &,
    gsl::span<int const> syms);
template uint32_t representative_subset<uint32_t, GroupActionLookup<uint32_t>>(
    uint32_t state, GroupActionLookup<uint32_t> const &,
    gsl::span<int const> syms);
template uint64_t representative_subset<uint64_t, GroupActionLookup<uint64_t>>(
    uint64_t state, GroupActionLookup<uint64_t> const &,
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
representative_sym<uint16_t, GroupAction>(uint16_t state, GroupAction const &);
template std::pair<uint32_t, int>
representative_sym<uint32_t, GroupAction>(uint32_t state, GroupAction const &);
template std::pair<uint64_t, int>
representative_sym<uint64_t, GroupAction>(uint64_t state, GroupAction const &);
template std::pair<uint16_t, int>
representative_sym<uint16_t, GroupActionLookup<uint16_t>>(
    uint16_t state, GroupActionLookup<uint16_t> const &);
template std::pair<uint32_t, int>
representative_sym<uint32_t, GroupActionLookup<uint32_t>>(
    uint32_t state, GroupActionLookup<uint32_t> const &);
template std::pair<uint64_t, int>
representative_sym<uint64_t, GroupActionLookup<uint64_t>>(
    uint64_t state, GroupActionLookup<uint64_t> const &);

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
representative_sym_subset<uint16_t, GroupAction>(uint16_t state,
                                                 GroupAction const &,
                                                 gsl::span<int const> syms);
template std::pair<uint32_t, int>
representative_sym_subset<uint32_t, GroupAction>(uint32_t state,
                                                 GroupAction const &,
                                                 gsl::span<int const> syms);
template std::pair<uint64_t, int>
representative_sym_subset<uint64_t, GroupAction>(uint64_t state,
                                                 GroupAction const &,
                                                 gsl::span<int const> syms);
template std::pair<uint16_t, int>
representative_sym_subset<uint16_t, GroupActionLookup<uint16_t>>(
    uint16_t state, GroupActionLookup<uint16_t> const &,
    gsl::span<int const> syms);
template std::pair<uint32_t, int>
representative_sym_subset<uint32_t, GroupActionLookup<uint32_t>>(
    uint32_t state, GroupActionLookup<uint32_t> const &,
    gsl::span<int const> syms);
template std::pair<uint64_t, int>
representative_sym_subset<uint64_t, GroupActionLookup<uint64_t>>(
    uint64_t state, GroupActionLookup<uint64_t> const &,
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
representative_syms<uint16_t, GroupAction>(uint16_t state, GroupAction const &);
template std::pair<uint32_t, std::vector<int>>
representative_syms<uint32_t, GroupAction>(uint32_t state, GroupAction const &);
template std::pair<uint64_t, std::vector<int>>
representative_syms<uint64_t, GroupAction>(uint64_t state, GroupAction const &);
template std::pair<uint16_t, std::vector<int>>
representative_syms<uint16_t, GroupActionLookup<uint16_t>>(
    uint16_t state, GroupActionLookup<uint16_t> const &);
template std::pair<uint32_t, std::vector<int>>
representative_syms<uint32_t, GroupActionLookup<uint32_t>>(
    uint32_t state, GroupActionLookup<uint32_t> const &);
template std::pair<uint64_t, std::vector<int>>
representative_syms<uint64_t, GroupActionLookup<uint64_t>>(
    uint64_t state, GroupActionLookup<uint64_t> const &);

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
mapping_syms<uint16_t, GroupAction>(uint16_t, uint16_t, GroupAction const &);
template std::vector<int>
mapping_syms<uint32_t, GroupAction>(uint32_t, uint32_t, GroupAction const &);
template std::vector<int>
mapping_syms<uint64_t, GroupAction>(uint64_t, uint64_t, GroupAction const &);
template std::vector<int> mapping_syms<uint16_t, GroupActionLookup<uint16_t>>(
    uint16_t, uint16_t, GroupActionLookup<uint16_t> const &);
template std::vector<int> mapping_syms<uint32_t, GroupActionLookup<uint32_t>>(
    uint32_t, uint32_t, GroupActionLookup<uint32_t> const &);
template std::vector<int> mapping_syms<uint64_t, GroupActionLookup<uint64_t>>(
    uint64_t, uint64_t, GroupActionLookup<uint64_t> const &);

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

template double norm<uint16_t, GroupAction>(uint16_t, GroupAction const &,
                                            Representation const &);
template double norm<uint32_t, GroupAction>(uint32_t, GroupAction const &,
                                            Representation const &);
template double norm<uint64_t, GroupAction>(uint64_t, GroupAction const &,
                                            Representation const &);
template double norm<uint16_t, GroupActionLookup<uint16_t>>(
    uint16_t, GroupActionLookup<uint16_t> const &, Representation const &);
template double norm<uint32_t, GroupActionLookup<uint32_t>>(
    uint32_t, GroupActionLookup<uint32_t> const &, Representation const &);
template double norm<uint64_t, GroupActionLookup<uint64_t>>(
    uint64_t, GroupActionLookup<uint64_t> const &, Representation const &);

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

template double norm_fermionic<uint16_t, GroupAction>(uint16_t state,
                                                      GroupAction const &,
                                                      Representation const &);
template double norm_fermionic<uint32_t, GroupAction>(uint32_t state,
                                                      GroupAction const &,
                                                      Representation const &);
template double norm_fermionic<uint64_t, GroupAction>(uint64_t state,
                                                      GroupAction const &,
                                                      Representation const &);
template double norm_fermionic<uint16_t, GroupActionLookup<uint16_t>>(
    uint16_t state, GroupActionLookup<uint16_t> const &,
    Representation const &);
template double norm_fermionic<uint32_t, GroupActionLookup<uint32_t>>(
    uint32_t state, GroupActionLookup<uint32_t> const &,
    Representation const &);
template double norm_fermionic<uint64_t, GroupActionLookup<uint64_t>>(
    uint64_t state, GroupActionLookup<uint64_t> const &,
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

template double norm_electron<uint16_t, GroupAction>(uint16_t, uint16_t,
                                                     GroupAction const &,
                                                     Representation const &);
template double norm_electron<uint32_t, GroupAction>(uint32_t, uint32_t,
                                                     GroupAction const &,
                                                     Representation const &);
template double norm_electron<uint64_t, GroupAction>(uint64_t, uint64_t,
                                                     GroupAction const &,
                                                     Representation const &);
template double norm_electron<uint16_t, GroupActionLookup<uint16_t>>(
    uint16_t, uint16_t, GroupActionLookup<uint16_t> const &,
    Representation const &);
template double norm_electron<uint32_t, GroupActionLookup<uint32_t>>(
    uint32_t, uint32_t, GroupActionLookup<uint32_t> const &,
    Representation const &);
template double norm_electron<uint64_t, GroupActionLookup<uint64_t>>(
    uint64_t, uint64_t, GroupActionLookup<uint64_t> const &,
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

template double norm_electron_subset<uint16_t, GroupAction>(
    uint16_t, uint16_t, GroupAction const &, Representation const &,
    gsl::span<int const>);
template double norm_electron_subset<uint32_t, GroupAction>(
    uint32_t, uint32_t, GroupAction const &, Representation const &,
    gsl::span<int const>);
template double norm_electron_subset<uint64_t, GroupAction>(
    uint64_t, uint64_t, GroupAction const &, Representation const &,
    gsl::span<int const>);
template double norm_electron_subset<uint16_t, GroupActionLookup<uint16_t>>(
    uint16_t, uint16_t, GroupActionLookup<uint16_t> const &,
    Representation const &, gsl::span<int const>);
template double norm_electron_subset<uint32_t, GroupActionLookup<uint32_t>>(
    uint32_t, uint32_t, GroupActionLookup<uint32_t> const &,
    Representation const &, gsl::span<int const>);
template double norm_electron_subset<uint64_t, GroupActionLookup<uint64_t>>(
    uint64_t, uint64_t, GroupActionLookup<uint64_t> const &,
    Representation const &, gsl::span<int const>);

} // namespace hydra::symmetries
