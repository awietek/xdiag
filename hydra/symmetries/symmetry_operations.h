#pragma once

#include <lila/external/gsl/span>

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

using span_size_t = gsl::span<int const>::size_type;

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
                            gsl::span<int const> syms);

// Computes the representative of state AND the symmetry that yields it
template <typename bit_t, class GroupAction>
std::pair<bit_t, int> representative_sym(bit_t state,
                                         GroupAction const &group_action);

// Computes the representative/symmetry using only a specified subset of
// symmetries
template <typename bit_t, class GroupAction>
std::pair<bit_t, int> representative_sym_subset(bit_t state,
                                                GroupAction const &group_action,
                                                gsl::span<int const> syms);

// Computes the symmetries, which map origin to target
template <typename bit_t, class GroupAction>
std::vector<int> mapping_syms(bit_t origin, bit_t target,
                              GroupAction const &group_action);

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
                     Representation const &irrep, gsl::span<int const> syms);

} // namespace hydra::symmetries
