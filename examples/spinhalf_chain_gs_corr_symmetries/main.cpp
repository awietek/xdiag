#include <xdiag/all.hpp>

using namespace xdiag;

std::complex<double> C_N_character(int N, int k, int p) {
  return std::exp(((2 * k * p) * M_PI / (double)N) *
                  std::complex<double>(0, 1.));
}

int main() try {
  // periodic N-site Heisenberg antiferromagnet
  int N = 14;
  auto H = OpSum();
  for (int i = 0; i < N; i++) {
    H += "J" * Op("SdotS", {i, (i + 1) % N});
  }
  H["J"] = 1.0;

  // define shift by one site
  std::vector<int64_t> T_perm(N, 0);
  for (int i = 0; i < N - 1; i++) {
    T_perm[i] = i + 1;
  }
  auto T = Permutation(T_perm);

  // define cyclic group
  std::vector<Permutation> perms(N);
  for (int i = 0; i < N; i++) {
    perms[i] = pow(T, i);
  }
  auto C_N = PermutationGroup(perms);

  // define irreps from character tables, labelled by momentum (2pi/N ×) k
  std::vector<Representation> irreps(N);
  std::vector<std::complex<double>> characters(N, 0.0);
  for (int k = 0; k < N; k++) {
    for (int p = 0; p < N; p++) {
      characters[p] = C_N_character(N, k, p);
    };
    irreps[k] = Representation(C_N, characters);
  }

  // check all translation symmetry blocks
  double e0 = 0.0;
  State psi0;
  for (int k = 0; k < N; k++) {
    auto block = Spinhalf(N, irreps[k]);
    auto [e0k, psik] = eig0(H, block);
    if (e0k < e0) {
      e0 = e0k;
      psi0 = psik;
    }
  }

  // we measure the correlation (which needs to besymmetrized)
  auto corr_op = Op("SdotS", {0, N / 2 - 1});
  auto corr = inner(symmetrize(corr_op, C_N), psi0);

  Log("Ground state energy: {:.12f}", e0);
  Log("Ground state correlator: {:.12f}", corr);
} catch (Error e) {
  error_trace(e);
}
