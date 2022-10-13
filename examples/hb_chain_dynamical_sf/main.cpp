#include <hydra/all.h>

int main() {
  using namespace hydra;
  using namespace arma;

  int n_sites = 16;
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

  // Create the irreps at momenta k
  std::vector<Representation> irreps;
  for (int k = 0; k < n_sites; ++k) {
    complex phase = exp(2i * pi * k / n_sites);
    auto irrep = generated_irrep(perm, phase);
    irreps.push_back(irrep);
  }

  // compute groundstate (known to be at k=0)
  Log("Computing ground state ...");
  auto block = Spinhalf(n_sites, n_up, group, irreps[0]);
  auto gs = groundstate(bonds, block);
  Log("done.");

  for (int q = 0; q < n_sites; ++q) {

    Log("Dynamical Lanczos iterations for q={}", q);

    // Compute S(q) |g.s.>
    auto S_of_q = symmetrized_operator(Bond("SZ", 0), group, irreps[q]);
    auto block_q = Spinhalf(n_sites, n_up, group, irreps[q]);
    auto v0 = zero_state(block_q);
    apply(S_of_q, gs, v0);

    double nrm = norm(v0);
    v0 /= nrm;

    // Perform 500 Lanczos iterations for dynamics starting from v0
    auto tmat = lanczos_eigenvalues_inplace(bonds, block_q, v0.vector(), 0, 0.,
                                            200, 1e-7);
    auto alphas = tmat.alphas();
    auto betas = tmat.betas();

    // Write alphas, betas, and norm to file for further processing
    std::stringstream sstr;
    sstr << ".N." << n_sites << ".nup." << n_up << ".q." << q << ".txt";
    alphas.save(std::string("outfiles/alphas") + sstr.str(), raw_ascii);
    betas.save(std::string("outfiles/betas") + sstr.str(), raw_ascii);
    std::ofstream of;
    of.open(std::string("outfiles/norm") + sstr.str());
    of << nrm;
    of.close();
  }

  return EXIT_SUCCESS;
}
