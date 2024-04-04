#pragma once
#include <vector>

#include <xdiag/extern/gsl/span>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/symmetries/operations/fermi_sign.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::symmetries {

// Computes stabilizing symmetries of "bits" when group is applied
template <typename bit_t, class GroupAction>
inline std::vector<int64_t> stabilizer_symmetries(bit_t bits,
                                                  GroupAction const &group) {
  std::vector<int64_t> stable_syms;
  for (int64_t sym = 0; sym < group.n_symmetries(); ++sym)
    if (group.apply(sym, bits) == bits)
      stable_syms.push_back(sym);
  return stable_syms;
}

// Computes the representative (smallest integer value) of "state"
// in orbit given by group_action
template <typename bit_t, class GroupAction>
inline bit_t representative(bit_t state, GroupAction const &group_action) {
  if (group_action.n_symmetries() == 0)
    return state;

  bit_t rep = std::numeric_limits<bit_t>::max();
  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < rep) {
      rep = tstate;
    }
  }
  return rep;
}

// determine whether a state is a representative
template <typename bit_t, class GroupAction>
inline bool is_representative(bit_t state, GroupAction const &group_action) {

  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < state) {
      return false;
    }
  }
  return true;
}

// Computes the representative using only a specified subset of symmetries
template <typename bit_t, class GroupAction>
inline bit_t representative_subset(bit_t state, GroupAction const &group_action,
                                   gsl::span<int64_t const> syms) {
  if (syms.size() == 0)
    return state;
  bit_t rep = std::numeric_limits<bit_t>::max();
  for (int64_t sym : syms) {
    assert(sym < group_action.n_symmetries());
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < rep) {
      rep = tstate;
    }
  }
  return rep;
}

// Computes the representative of state AND the symmetry that yields it
template <typename bit_t, class GroupAction>
inline std::pair<bit_t, int64_t>
representative_sym(bit_t state, GroupAction const &group_action) {
  if (group_action.n_symmetries() == 0)
    return {state, 0};
  bit_t rep = std::numeric_limits<bit_t>::max();
  int64_t rep_sym = 0;
  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < rep) {
      rep = tstate;
      rep_sym = sym;
    }
  }
  return {rep, rep_sym};
}

// Computes the representative/symmetry using only a specified subset of
// symmetries
template <typename bit_t, class GroupAction>
inline std::pair<bit_t, int64_t>
representative_sym_subset(bit_t state, GroupAction const &group_action,
                          gsl::span<int64_t const> syms) {
  if (syms.size() == 0)
    return {state, 0};
  bit_t rep = std::numeric_limits<bit_t>::max();
  int64_t rep_sym = 0;
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

// Computes the representative of state and all symmetries that yield it
template <typename bit_t, class GroupAction>
inline std::pair<bit_t, std::vector<int64_t>>
representative_syms(bit_t state, GroupAction const &group_action) {
  if (group_action.n_symmetries() == 0)
    return {state, std::vector<int64_t>()};

  bit_t rep = representative(state, group_action);
  std::vector<int64_t> rep_syms;
  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate == rep) {
      rep_syms.push_back(sym);
    }
  }
  return {rep, rep_syms};
}

// Computes the symmetries, which map origin to target
template <typename bit_t, class GroupAction>
inline std::vector<int64_t> mapping_syms(bit_t origin, bit_t target,
                                         GroupAction const &group_action) {
  std::vector<int64_t> syms;
  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, origin);
    if (tstate == target) {
      syms.push_back(sym);
    }
  }
  return syms;
}

// Computes the norm of a symmetrized state
template <typename bit_t, class GroupAction>
double norm(bit_t state, GroupAction const &group_action,
            Representation const &irrep) {
  complex amplitude = 0.0;
  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate == state) {
      amplitude += irrep.character(sym);
    }
  }
  return std::sqrt(std::abs(amplitude));
}

// Computes the norm of a symmetrized state with fermions
template <typename bit_t, class GroupAction>
inline double norm_fermionic(bit_t state, GroupAction const &group_action,
                             Representation const &irrep) {
  complex amplitude = 0.0;
  int64_t n_sites = group_action.n_sites();
  auto work = fermi_work(n_sites);
  auto const &group = group_action.permutation_group();

  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {

    auto const &perm = group[sym];

    bit_t tstate = group_action.apply(sym, state);
    if (tstate == state) {
      if (fermi_bool_of_permutation(state, perm, work)) {
        amplitude -= irrep.character(sym);
      } else {
        amplitude += irrep.character(sym);
      }
    }
  }
  return std::sqrt(std::abs(amplitude));
}

// Computes the norm of a symmetrized state with up/dn electrons
template <typename bit_t, class GroupAction>
inline double norm_electron(bit_t ups, bit_t dns,
                            GroupAction const &group_action,
                            Representation const &irrep) {
  assert(group_action.n_symmetries() == irrep.size());
  complex amplitude = 0.0;
  int64_t n_sites = group_action.n_sites();
  auto work = fermi_work(n_sites);
  auto const &group = group_action.permutation_group();

  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {

    auto const &perm = group[sym];
    bit_t tups = group_action.apply(sym, ups);

    if (tups == ups) {
      bool fermi_bool_ups = fermi_bool_of_permutation(ups, perm, work);

      bit_t tdns = group_action.apply(sym, dns);

      if (tdns == dns) {
        bool fermi_bool_dns = fermi_bool_of_permutation(dns, perm, work);

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

// Computes the norm of a symmetrized state with up/dn electrons (subset of
// syms)
template <typename bit_t, class GroupAction>
inline double norm_electron_subset(bit_t ups, bit_t dns,
                                   GroupAction const &group_action,
                                   Representation const &irrep,
                                   gsl::span<int64_t const> syms) {
  assert(group_action.n_symmetries() == irrep.size());
  complex amplitude = 0.0;
  int64_t n_sites = group_action.n_sites();
  auto work = fermi_work(n_sites);
  auto const &group = group_action.permutation_group();

  for (int64_t sym : syms) {
    assert(sym < group_action.n_symmetries());

    auto const &perm = group[sym];

    bit_t tups = group_action.apply(sym, ups);

    if (tups == ups) {
      bool fermi_bool_ups = fermi_bool_of_permutation(ups, perm, work);

      bit_t tdns = group_action.apply(sym, dns);

      if (tdns == dns) {
        bool fermi_bool_dns = fermi_bool_of_permutation(dns, perm, work);
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

} // namespace xdiag::symmetries
