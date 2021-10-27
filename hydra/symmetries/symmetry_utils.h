#pragma once

#include <hydra/combinatorics/combinations.h>

#include <hydra/symmetries/fermi_sign.h>
#include <hydra/symmetries/representation.h>

#include <string>
#include <utility>
#include <vector>

namespace hydra {
std::vector<std::vector<int>> read_permutations(std::string filename);
}

namespace hydra::utils {

bool is_valid_permutation(int n_sites, const int *permutation);

template <class bit_t>
bit_t apply_permutation(bit_t state, int n_sites, const int *permutation);

template <class bit_t, class GroupAction>
std::vector<int> stabilizer_symmetries(bit_t bits, GroupAction &&group) {
  std::vector<int> stable_syms;
  for (int sym = 0; sym < group.n_symmetries(); ++sym)
    if (group.apply(sym, bits) == bits)
      stable_syms.push_back(sym);
  return stable_syms;
}

template <class bit_t, class GroupAction>
bit_t representative(bit_t state, GroupAction &&group_action) {
  bit_t rep = std::numeric_limits<bit_t>::max();
  for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate < rep) {
      rep = tstate;
    }
  }
  return rep;
}

template <class bit_t, class GroupAction>
std::pair<bit_t, int> representative_sym(bit_t state,
                                         GroupAction &&group_action) {
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

template <class bit_t, class GroupAction>
double compute_norm(bit_t state, GroupAction &&group_action,
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

template <class bit_t, class GroupAction>
double compute_norm_fermionic(bit_t state, GroupAction &&group_action,
                              Representation const &irrep, int *work) {
  // "work" needs to be allocated of size n_sites

  complex amplitude = 0.0;
  int n_sites = group_action.n_sites();
  const int *sym_ptr = group_action.permutation_array().data();
  for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate == state) {
      amplitude += irrep.character(sym) *
                   fermi_sign_of_permutation(state, sym_ptr, work);
    }
    sym_ptr += n_sites;
  }
  return std::sqrt(std::abs(amplitude));
}

template <class bit_t, class GroupAction, class LinTable>
void fill_reps_idces_syms_limits(
    int n_sites, int n_par, GroupAction &&group_action, LinTable &&lintable,
    std::vector<bit_t> &reps, std::vector<idx_t> &idces, std::vector<int> &syms,
    std::vector<std::pair<idx_t, idx_t>> &sym_limits) {

  using combinatorics::Combinations;

  // Compute all representatives
  idx_t idx = 0;
  for (bit_t state : Combinations<bit_t>(n_sites, n_par)) {
    bit_t rep = utils::representative(state, group_action);
    if (rep == state) {
      idces[idx] = reps.size();
      reps.push_back(rep);
    }
    ++idx;
  }

  // Compute indices of up-representatives and stabilizer symmetries
  idx = 0;
  for (bit_t state : Combinations<bit_t>(n_sites, n_par)) {
    bit_t rep = utils::representative(state, group_action);
    idces[idx] = idces[lintable.index(rep)];

    // Determine the symmetries that yield the up-representative
    idx_t begin = syms.size();
    for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
      if (group_action.apply(sym, state) == rep)
        syms.push_back(sym);
    }
    idx_t end = syms.size();
    sym_limits[idx] = {begin, end};
    ++idx;
  }
}

template <class bit_t, class GroupAction>
void fill_fermi_bool_table(GroupAction &&group_action, int npar,
                           std::vector<bool> &fermi_bool_table) {
  using combinatorics::Combinations;

  int n_sites = group_action.n_sites();
  int n_symmetries = group_action.n_symmetries();
  size_t raw_size = combinatorics::binomial(n_sites, npar);
  std::vector<int> fermi_work(n_sites, 0);
  const int *sym_ptr = group_action.permutation_array().data();
  for (int sym = 0; sym < n_symmetries; ++sym) {
    idx_t idx = 0;
    for (bit_t state : Combinations<bit_t>(n_sites, npar)) {
      fermi_bool_table[sym * raw_size + idx] =
          utils::fermi_bool_of_permutation(state, sym_ptr, fermi_work.data());
      ++idx;
    }
    sym_ptr += n_sites;
  }
}

} // namespace hydra::utils
