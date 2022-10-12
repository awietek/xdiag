#include <hydra/all.h>

int main() {
  using namespace hydra;

  Log.set_verbosity(1);

  int n_sites = 24;
  int n_up = n_sites / 2;

  // Create nearest-neighbor Heisenberg model
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

  // Determine the ground state block
  double e0 = 0.;
  int k0 = 0;
  for (int k = 0; k < n_sites; ++k) {

    // Create the irrep with momentum "k"
    complex phase = exp(2i * pi * k / n_sites);
    auto irrep = generated_irrep(perm, phase);

    // compute the ground state energy in this sector
    auto block = Spinhalf(n_sites, n_up, group, irrep);
    auto e0_k = eig0(bonds, block);

    Log("k: {}, e0: {}", k, e0_k);

    if (e0_k < e0) {
      e0 = e0_k;
      k0 = k;
    }
  }

  // Create the irrep with momentum "k0"
  complex phase = exp(2i * pi * k0 / n_sites);
  auto irrep = generated_irrep(perm, phase);

  // Compute ground state
  auto block = Spinhalf(n_sites, n_up, group, irrep);
  auto gs = groundstate(bonds, block);
  HydraPrint(gs);

  return EXIT_SUCCESS;
}
