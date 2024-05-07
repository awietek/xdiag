#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {

  int n_sites = 16;
  int nup = n_sites / 2;

  // Define the nearest-neighbor Heisenberg Hamiltonian
  BondList bonds;
  for (int i = 0; i < n_sites; ++i) {
    bonds << Bond("HB", "J", {i, (i + 1) % n_sites});
  }
  bonds["J"] = 1.0;

  // Create the permutation group
  std::vector<int> translation;
  for (int s = 0; s < n_sites; ++s) {
    translation.push_back((s + 1) % n_sites);
  }
  Permutation perm(translation);
  auto group = generated_group(perm);

  // Loop over momenta k
  for (int k = 0; k < n_sites; ++k) {

    // Create irrep with momentum k
    complex phase = exp(2i * pi * k / n_sites);
    auto irrep = generated_irrep(perm, phase);

    // Compute the groundstate energy
    auto block = Spinhalf(n_sites, nup, group, irrep);
    double e0 = eigval0(bonds, block);
    Log("k: {}, e0: {}", k, e0);
  }
} catch (Error e) {
  error_trace(e);
}
