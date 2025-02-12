#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {

  int nsites = 16;
  int nup = nsites / 2;

  // Define the nearest-neighbor Heisenberg Hamiltonian
  OpSum ops;
  for (int i = 0; i < nsites; ++i) {
    ops += Op("HB", "J", {i, (i + 1) % nsites});
  }
  ops["J"] = 1.0;

  // Create the permutation group
  std::vector<int> translation;
  for (int s = 0; s < nsites; ++s) {
    translation.push_back((s + 1) % nsites);
  }
  Permutation perm(translation);
  auto group = generated_group(perm);

  // Loop over momenta k
  for (int k = 0; k < nsites; ++k) {

    // Create irrep with momentum k
    complex phase = exp(2i * pi * k / nsites);
    auto irrep = generated_irrep(perm, phase);

    // Compute the groundstate energy
    auto block = Spinhalf(nsites, nup, group, irrep);
    double e0 = eigval0(ops, block);
    Log("k: {}, e0: {}", k, e0);
  }
} catch (Error e) {
  error_trace(e);
}
