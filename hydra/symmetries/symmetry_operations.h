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

namespace hydra::symmetries {

bool is_valid_permutation(int n_sites, const int *permutation);

template <typename bit_t>
bit_t apply_permutation(bit_t state, int n_sites, const int *permutation);

// Computes stabilizing symmetries of "bits" when group is applied
template <typename bit_t, class GroupAction>
std::vector<int> stabilizer_symmetries(bit_t bits, GroupAction const &group);

// Computes the representative (smallest integer value) of "state"
// in orbit given by group_action
template <typename bit_t, class GroupAction>
bit_t representative(bit_t state, GroupAction const &group_action);

// Computes the representative using only a specified subset of symmetries
template <typename bit_t, class GroupAction>
bit_t representative_subset(bit_t state, GroupAction const &group_action,
                            std::vector<int> const &syms);

// Computes the representative of state AND the symmetry that yields it
template <typename bit_t, class GroupAction>
std::pair<bit_t, int> representative_sym(bit_t state,
                                         GroupAction const &group_action);

// Computes the representative/symmetry using only a specified subset of
// symmetries
template <typename bit_t, class GroupAction>
std::pair<bit_t, int> representative_sym_subset(bit_t state,
                                                GroupAction const &group_action,
                                                std::vector<int> const &syms);

// Computes the norm of a symmetrized state
template <typename bit_t, class GroupAction>
double norm(bit_t state, GroupAction const &group_action,
            Representation const &irrep);

// Computes the norm of a symmetrized state with fermions
template <typename bit_t, class GroupAction>
double norm_fermionic(bit_t state, GroupAction const &group_action,
                      Representation const &irrep);

// Computes the norm of a symmetrized state with up/dn electrons
template <typename bit_t, class GroupAction>
double norm_electron(bit_t ups, bit_t dns, GroupAction const &group_action,
                     Representation const &irrep);

// Computes the norm of a symmetrized state with up/dn electrons (subset of
// syms)
template <typename bit_t, class GroupAction>
double
norm_electron_subset(bit_t ups, bit_t dns, GroupAction const &group_action,
                     Representation const &irrep, std::vector<int> const &syms);

// Create a lists of representatives and their symmetries
template <typename bit_t, class GroupAction, class LinTable>
std::tuple<std::vector<bit_t>, std::vector<idx_t>, std::vector<int>,
           std::vector<std::pair<idx_t, idx_t>>>
representatives_indices_symmetries_limits(int npar,
                                          GroupAction const &group_action,
                                          LinTable const &lintable);

// Creates the table of Fermi signs for a given particle number
template <typename bit_t, class GroupAction>
std::vector<bool> fermi_bool_table(int npar, GroupAction const &group_action);

} // namespace hydra::symmetries
