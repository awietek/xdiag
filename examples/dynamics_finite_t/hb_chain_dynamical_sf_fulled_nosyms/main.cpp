#include <filesystem>
#include <xdiag/all.hpp>

int main() {
  using namespace xdiag;
  using namespace arma;
  using fmt::format;
  using hdf5_opts::append;

  say_hello();
  
  // Parse input arguments
  int nsites = 12;

  // Define directory / file to store output data
  std::string outdir = format("outfiles/N.{}", nsites);
  std::string outfile = format("{}/outfile.N.{}.h5", outdir, nsites);
  std::filesystem::create_directories(outdir);

  // Create nearest-neighbor Heisenberg model
  OpSum ops;
  for (int s = 0; s < nsites; ++s) {
    ops << Op("HB", "J", {s, (s + 1) % nsites});
  }
  ops["J"] = 1.0;

  auto block = Spinhalf(nsites);

  // Compute eigendecomposition of Hamiltonian
  Log("Creating H");
  mat H = matrix(ops, block);

  Log("Diagonalizing H");
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, H);
  eigval.save(hdf5_name(outfile, "eigenvalues", append));

  // Loop over different momenta q
  for (int q = 0; q <= nsites; ++q) {

    Log("Computing <n|S(q)|m> (q={})", q);
    
    // Create S(q) operator
    OpSum S_of_q_ops;
    for (int s = 0; s < nsites; ++s) {
      complex phase = exp(2i * pi * q * s / nsites);
      S_of_q_ops << Op("Sz", phase / nsites, s);
    }

    // Compute matrix elements of S(q)
    cx_mat S_of_q = matrixC(S_of_q_ops, block);
    cx_mat S_of_q_eig = eigvec.t() * S_of_q * eigvec;
    S_of_q_eig.save(hdf5_name(outfile, format("S_of_q_{}_eig", q), append));
  }

  return EXIT_SUCCESS;
}
