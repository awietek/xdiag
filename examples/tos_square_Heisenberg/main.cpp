#include <xdiag/all.hpp>

int main(int argc, char **argv) {
  using namespace xdiag;
  using namespace arma;
  using fmt::format;

  say_hello();

  // Parse input arguments
  assert(argc == 6);
  int n_sites = atoi(argv[1]);              // number of sites
  int n_up = atoi(argv[2]);                 // number of upspins
  std::string kname = std::string(argv[3]); // momentum k
  double J1 = atof(argv[4]);
  int seed = atoi(argv[5]);

  Log("Diagonalizing H in block nup: {}, k: {}", n_up, kname);

  auto lfile = FileToml(format("square.{}.J1.fsl.pbc.toml", n_sites));
  std::string ofilename =
      format("outfile.square.{}.J1.{:.2f}.nup.{}.k.{}.seed.{}.h5", n_sites, J1,
             n_up, kname, seed);
  auto ofile = FileH5(ofilename, "w!");

  OpSum ops = read_opsum(lfile, "Interactions");
  ops["J1"] = J1;
  auto irrep = read_representation(lfile, kname);

  Log("Creating block ...");
  tic();
  auto block = Spinhalf(n_sites, n_up, irrep);
  toc();
  Log("Dimension: {}", block.size());

  Log("Running Lanczos ...");
  tic();
  int n_eig_to_converge = 2;
  int max_iterations = 40;
  auto tmat = eigvals_lanczos(ops, block, n_eig_to_converge, 1e-12,
                              max_iterations, 1e-7, seed);
  toc();

  ofile["Alphas"] = tmat.alphas;
  ofile["Betas"] = tmat.betas;
  ofile["Eigenvalues"] = tmat.eigenvalues;
  ofile["Dimension"] = block.size();

  return EXIT_SUCCESS;
}
