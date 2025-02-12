#include <xdiag/all.hpp>

int main(int argc, char **argv) {
  using namespace xdiag;
  using namespace arma;
  using fmt::format;

  say_hello();
  
  // Parse input arguments
  assert(argc == 6);
  int nsites = atoi(argv[1]);              // number of sites
  int nup = atoi(argv[2]);                 // number of upspins
  std::string kname = std::string(argv[3]); // momentum k
  double J = atof(argv[4]);
  double Jd = atof(argv[5]);

  Log("Diagonalizing H in block nup: {}, k: {}", nup, kname);

  auto lfile = FileToml(format("checkerboard.{}.JJD.toml", nsites), 'r');
  auto ofile = FileH5(
      format(
          "outfiles/outfile.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.nup.{}.k.{}.h5",
          nsites, J, Jd, nup, kname),
      "w!");

  auto ops = OpSum(lfile["Interactions"]);
  ops["J"] = J;
  ops["JD"] = Jd;
  auto group = PermutationGroup(lfile["Symmetries"]);
  auto irrep_k = Representation(lfile[kname]);

  Log("Creating block ...");
  tic();
  auto block_k = Spinhalf(nsites, nup, group, irrep_k);
  toc();
  Log("Dimension: {}", block_k.size());

  Log("Creating matrix ...");
  tic();
  auto H = matrix(ops, block_k);
  toc();

  Log("Diagonalizing ...");
  tic();
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, H);
  ofile["Eigenvalues"] = eigval;
  toc();

  return EXIT_SUCCESS;
}
