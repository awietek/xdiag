#include <xdiag/all.hpp>

int main() {
  using namespace xdiag;
  using namespace arma;

  say_hello();
  
  int nsites = 12;
  int nup = nsites / 2;

  // Create nearest-neighbor Heisenberg model
  OpSum ops;
  for (int s = 0; s < nsites; ++s) {
    ops << Op("HB", "J", {s, (s + 1) % nsites});
  }
  ops["J"] = 1.0;

  // compute groundstate (known to be at k=0)
  Log("Computing ground state ...");
  auto block = Spinhalf(nsites, nup);
  // auto gs = groundstate(ops, block);
  auto resgs = eigs_lanczos(ops, block);
  auto gs = resgs.eigenvectors.col(0);
  gs.make_complex();
  Log("done.");

  for (int q = 0; q < nsites; ++q) {

    Log("Dynamical Lanczos iterations for q={}", q);

    complex phase = exp(2i * pi * q / nsites);
    OpSum S_of_q;
    for (int s = 0; s < nsites; ++s) {
      S_of_q << Op("Sz", pow(phase, s) / nsites, s);
    }

    // Compute S(q) |g.s.>
    auto v0 = State(block);
    v0.make_complex();
    apply(S_of_q, gs, v0);
    double nrm = norm(v0);
    v0 /= nrm;

    // Perform 200 Lanczos iterations for dynamics starting from v0
    auto res = eigvals_lanczos(ops, block, v0, 1, 0., 200, true, 1e-7);
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
