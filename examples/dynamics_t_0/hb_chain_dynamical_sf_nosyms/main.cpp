#include <xdiag/all.hpp>

int main() {
  using namespace xdiag;
  using namespace arma;

  say_hello();
  
  int n_sites = 12;
  int n_up = n_sites / 2;

  // Create nearest-neighbor Heisenberg model
  OpSum ops;
  for (int s = 0; s < n_sites; ++s) {
    ops << Op("HB", "J", {s, (s + 1) % n_sites});
  }
  ops["J"] = 1.0;

  // compute groundstate (known to be at k=0)
  Log("Computing ground state ...");
  auto block = Spinhalf(n_sites, n_up);
  // auto gs = groundstate(ops, block);
  auto resgs = eigs_lanczos(ops, block);
  auto gs = resgs.eigenvectors.col(0);
  gs.make_complex();
  Log("done.");

  for (int q = 0; q < n_sites; ++q) {

    Log("Dynamical Lanczos iterations for q={}", q);

    complex phase = exp(2i * pi * q / n_sites);
    OpSum S_of_q;
    for (int s = 0; s < n_sites; ++s) {
      S_of_q << Op("Sz", pow(phase, s) / n_sites, s);
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
        fmt::format("outfiles/alphas.N.{}.nup.{}.q.{}.txt", n_sites, n_up, q),
        raw_ascii);
    betas.save(
        fmt::format("outfiles/betas.N.{}.nup.{}.q.{}.txt", n_sites, n_up, q),
        raw_ascii);
    std::ofstream of;
    of.open(
        fmt::format("outfiles/norm.N.{}.nup.{}.q.{}.txt", n_sites, n_up, q));
    of << nrm;
    of.close();
  }

  return EXIT_SUCCESS;
}
