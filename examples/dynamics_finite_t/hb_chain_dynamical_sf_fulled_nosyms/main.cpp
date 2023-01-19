#include <filesystem>
#include <hydra/all.h>

int main() {
  using namespace hydra;
  using namespace arma;
  using fmt::format;
  using hdf5_opts::append;

  // Parse input arguments
  int n_sites = 12;

  // Define directory / file to store output data
  std::string outdir = format("outfiles/N.{}", n_sites);
  std::string outfile = format("{}/outfile.N.{}.h5", outdir, n_sites);
  std::filesystem::create_directories(outdir);

  // Create nearest-neighbor Heisenberg model
  BondList bonds;
  for (int s = 0; s < n_sites; ++s) {
    bonds << Bond("HB", "J", {s, (s + 1) % n_sites});
  }
  bonds["J"] = 1.0;

  auto block = Spinhalf(n_sites);

  // Compute eigendecomposition of Hamiltonian
  auto H = matrix_real(block, bonds);
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, H);
  eigval.save(hdf5_name(outfile, format("eigenvalues", q), append));

  // Loop over different momenta q
  for (int q = 0; q <= n_sites; ++q) {

    // Create S(q) operator
    BondList S_of_q_bonds;
    for (int s = 0; s < n_sites; ++s) {
      complex phase = exp(2i * pi * q * s / n_sites);
      S_of_q_bonds << Bond("SZ", phase, s);
    }

    // Compute matrix elements of S(q)
    auto S_of_q = matrix_cplx(block, S_of_q_bonds);
    auto S_of_q_eig = eigvec.t() * S_of_q * eigvec;
    S_of_q_eig.save(hdf5_name(outfile, format("S_of_q_{}_eig", q), append));
  }

  return EXIT_SUCCESS;
}
