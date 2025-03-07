#include <xdiag/all.hpp>

using namespace xdiag;
using namespace std::complex_literals;
using fmt::format;

constexpr double XDIAG_PI = 3.14159265358979323846;

int main() try {
  say_hello();

  // IO
  std::string filename = XDIAG_DIRECTORY
      "/misc/data/examples_output/spinhalf_chain_structure_factor.h5";
  auto outfile = FileH5(filename, "w!");

  // Define the nearest-neighbor Heisenberg model
  int N = 16;
  int nup = N / 2;
  OpSum ops;
  for (int i = 0; i < N; ++i) {
    ops += "J" * Op("SdotS", {i, (i + 1) % N});
  }
  ops["J"] = 1.0;

  // Create the permutation group
  std::vector<int> translation;
  for (int s = 0; s < N; ++s) {
    translation.push_back((s + 1) % N);
  }
  Permutation perm(translation);
  std::vector<Permutation> perms;
  for (int s = 0; s < N; ++s) {
    perms.push_back(pow(perm, s));
  }
  auto group = PermutationGroup(perms);

  // Create the irreps at momenta k
  std::vector<Representation> irreps;
  for (int k = 0; k < N; ++k) {
    complex phase = exp(2i * XDIAG_PI * (double)k / (double)N);
    std::vector<complex> characters;
    for (int s = 0; s < N; ++s) {
      characters.push_back(pow(phase, s));
    }
    auto irrep = Representation(group, characters);
    irreps.push_back(irrep);
  }

  // compute groundstate (known to be at k=0)
  Log("Computing ground state ...");
  auto block = Spinhalf(N, nup, irreps[0]);
  auto [e0, gs] = eig0(ops, block);
  gs.make_complex();
  Log("done.");
  Log("Ground state energy: {:.12f}", e0);
  outfile["e0"] = e0;

  // loop through momenta
  for (int k = 0; k < N; ++k) {
    Log("S^zz(k,w) at momentum {}", k);
    auto S_q = symmetrize(Op("Sz", 0), irreps[k]);
    auto Av = apply(S_q, gs);
    auto nrm = norm(Av);
    Av /= nrm;
    XDIAG_SHOW(nrm);

    auto res = eigvals_lanczos_inplace(ops, Av);
    outfile[format("{}/norm", k)] = nrm;
    outfile[format("{}/alphas", k)] = res.alphas;
    outfile[format("{}/betas", k)] = res.betas;
    outfile[format("{}/eigs", k)] = res.eigenvalues;
  }

} catch (Error e) {
  error_trace(e);
}
