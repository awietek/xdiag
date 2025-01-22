#include <xdiag/all.hpp>

int main() {
  using namespace xdiag;
  using namespace arma;

  say_hello();
  
  int nsites = 16;
  int nup = nsites / 2;

  // Create nearest-neighbor Heisenberg model
  OpSum ops;
  for (int s = 0; s < nsites; ++s) {
    ops << Op("HB", "J", {s, (s + 1) % nsites});
  }
  ops["J"] = 1.0;

  // Create the permutation group
  std::vector<int> translation;
  for (int s = 0; s < nsites; ++s) {
    translation.push_back((s + 1) % nsites);
  }
  Permutation perm(translation);
  auto group = generated_group(perm);

  // Create the irreps at momenta k
  std::vector<Representation> irreps;
  for (int k = 0; k < nsites; ++k) {
    complex phase = exp(2i * pi * k / nsites);
    auto irrep = generated_irrep(perm, phase);
    irreps.push_back(irrep);
  }

  // compute groundstate (known to be at k=0)
  Log("Computing ground state ...");
  auto block = Spinhalf(nsites, nup, group, irreps[0]);
  auto [e0, gs] = eig0(ops, block);
  // auto gs = resgs.eigenvectors.col(0);
  gs.make_complex();
  Log("done.");

  for (int q = 0; q < nsites; ++q) {

    Log("Dynamical Lanczos iterations for q={}", q);

    // Compute S(q) |g.s.>
    auto S_of_q = symmetrized_operator(Op("Sz", 0), group, irreps[q]);
    auto block_q = Spinhalf(nsites, nup, group, irreps[q]);
    auto v0 = State(block_q);
    v0.make_complex();
    apply(S_of_q, gs, v0);


    double nrm = norm(v0);
    v0 /= nrm;

    // Perform 200 Lanczos iterations for dynamics starting from v0
    auto res = eigvals_lanczos(ops, block_q, v0, 1, 0., 200, true, 1e-7);
    auto alphas = res.alphas;
    auto betas = res.betas;

    // Write alphas, betas, and norm to file for further processing
    alphas.save(
        fmt::format("outfiles/alphas.N.{}.nup.{}.q.{}.txt", nsites, nup, q),
        raw_ascii);
    betas.save(
        fmt::format("outfiles/betas.N.{}.nup.{}.q.{}.txt", nsites, nup, q),
        raw_ascii);
    std::ofstream of;
    of.open(
        fmt::format("outfiles/norm.N.{}.nup.{}.q.{}.txt", nsites, nup, q));
    of << nrm;
    of.close();
  }

  return EXIT_SUCCESS;
}
