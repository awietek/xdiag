#include "generated_group.hpp"

#include <xdiag/symmetries/generated_group.hpp>

namespace xdiag {

PermutationGroup generated_group(Permutation const &generator) {
  return generated_group(std::vector<Permutation>{generator});
}
PermutationGroup generated_group(std::vector<Permutation> const &generators) {
  auto [group, irrep] = generated_group_irrep(
      generators, std::vector<complex>(generators.size(), 1.0));
  (void)irrep;
  return group;
}

Representation generated_irrep(Permutation const &generator, complex phase) {
  return generated_irrep(std::vector<Permutation>{generator}, {phase});
}
Representation generated_irrep(std::vector<Permutation> const &generators,
                               std::vector<complex> const &phases) {
  auto [group, irrep] = generated_group_irrep(generators, phases);
  (void)group;
  return irrep;
}

std::pair<PermutationGroup, Representation>
generated_group_irrep(std::vector<Permutation> const &generators,
                      std::vector<complex> const &phases) {
  if (generators.size() == 0) {
    return {PermutationGroup(), Representation()};
  }
  if (generators.size() != phases.size()) {
    Log.err("Error in generated_group_irrep: number of generators not equal to "
            "the number of phases");
  }

  int n_sites = generators[0].size();

  // Check whether all generators have the same number of sites
  for (std::size_t i = 0; i < generators.size(); ++i) {
    if (generators[i].size() != n_sites) {
      Log.err(
          "Error in generated_group: not all generators have the same size!");
    }
  }

  std::vector<Permutation> perms = {identity_permutation(n_sites)};
  std::vector<complex> characters = {1.0};

  bool nothing_to_add = false;
  while (!nothing_to_add) {

    auto updated_perms = perms;
    auto updated_characters = characters;

    bool new_perm_created = false;

    for (std::size_t p_idx = 0; p_idx < perms.size(); ++p_idx) {
      auto p = perms[p_idx];
      auto p_phase = characters[p_idx];

      for (std::size_t g_idx = 0; g_idx < generators.size(); ++g_idx) {
        auto g = generators[g_idx];
        auto g_phase = phases[g_idx];

        auto p_times_g = p * g;
        auto p_times_g_phase = p_phase * g_phase;
        if (std::find(updated_perms.begin(), updated_perms.end(), p_times_g) ==
            updated_perms.end()) {
          updated_perms.push_back(p_times_g);
          updated_characters.push_back(p_times_g_phase);
          new_perm_created = true;
        }
      }
    }
    if (!new_perm_created) {
      nothing_to_add = true;
    }
    perms = updated_perms;
    characters = updated_characters;
  }

  return {PermutationGroup(perms), Representation(characters)};
}

} // namespace xdiag
