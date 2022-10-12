#include <hydra/all.h>

int main() {
  using namespace hydra;
  
  int n_sites = 16;

  // Create nearest neighbor Heisenberg model
  BondList bonds;
  for (int s = 0; s < n_sites; ++s) {
    bonds << Bond("HB", "J", {s, (s + 1) % n_sites});
  }
  bonds["J"] = 1.0;

  // Create the permutation group
  std::vector<int> translation;
  for (int s = 0; s < n_sites; ++s) {
    translation.push_back((s + 1) % n_sites);
  }
  Permutation perm(translation);
  auto group = generated_group(perm);
  // HydraPrint(group);

  return EXIT_SUCCESS;
}
