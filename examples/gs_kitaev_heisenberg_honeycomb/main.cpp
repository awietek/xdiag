#include <filesystem>
#include <xdiag/all.hpp>

int main(int argc, char **argv) {
  using namespace xdiag;
  using namespace arma;
  using fmt::format;

  say_hello();

  // Parse input arguments
  assert(argc == 8);
  int n_sites = atoi(argv[1]);              // number of sites
  std::string kname = std::string(argv[2]); // momentum k
  double J = atof(argv[3]);
  double KX = atof(argv[4]);
  double KY = atof(argv[5]);
  double KZ = atof(argv[6]);
  int seed = atoi(argv[7]);

  Log("Diagonalizing H in block k: {}", kname);

  // auto lfile = FileToml(format("honeycomb.{}.HeisenbergKitaevGamma.fsl.toml", n_sites));
  auto lfile = FileToml(format("kitaev.{}.toml", n_sites));
  std::string odir = format("outfiles/seed.{}", seed);
  std::string ofilename = format(
      "outfile.honeycomb.{}.J.{:.2f}.KX.{:.2f}.KY.{:.2f}.KZ.{:.2f}.k.{}.seed.{}.h5",
      n_sites, J, KX, KY, KZ, kname, seed);
  std::filesystem::create_directories(odir);
  auto ofile = FileH5(format("{}/{}", odir, ofilename), "w!");

  // xdiag::OpSum ops = read_opsum(lfile, "Interactions");
  // ops["J"] = J;
  // ops["KX"] = KX;
  // ops["KY"] = KY;
  // ops["KZ"] = KZ;
  auto ops_read = read_opsum(lfile, "Interactions");
  cx_mat sx(mat({{0., 0.5}, {0.5, 0.}}), mat({{0., 0.}, {0., 0.}}));
  cx_mat sy(mat({{0., 0.}, {0., 0.}}), mat({{0., -0.5}, {0.5, 0.}}));
  cx_mat sz(mat({{0.5, 0.0}, {0.0, -0.5}}), mat({{0., 0.}, {0., 0.0}}));

  cx_mat sxsx = kron(sx, sx);
  cx_mat sysy = kron(sy, sy);
  cx_mat szsz = kron(sz, sz);
  cx_mat gsx = kron(sy, sz) + kron(sz, sy);
  cx_mat gsy = kron(sx, sz) + kron(sz, sx);
  cx_mat gsz = kron(sx, sy) + kron(sy, sx);

  auto ops = OpSum();
  for (auto [cpl, op] : ops_read) {
    std::string type = op.type();
    auto sites = op.sites();
    if (type == "KITAEVX") {
      ops += KX * Op("Matrix", sites, sxsx);
    } else if (type == "KITAEVY") {
      ops += KY * Op("Matrix", sites, sysy);
    } else if (type == "KITAEVZ") {
      ops += KZ * Op("Matrix", sites, szsz);
    } else if (type == "SdotS") {
      ops += J * Op("Matrix", sites, sxsx);
      ops += J * Op("Matrix", sites, sysy);
      ops += J * Op("Matrix", sites, szsz);
    }
  }
  auto irrep = read_representation(lfile, kname);

  Log("Creating block ...");
  tic();
  auto block = Spinhalf(n_sites, irrep);
  toc();
  Log("Dimension: {}", block.size());

  Log("Running Lanczos ...");
  tic();
  int n_eig_to_converge = 2;
  int max_iterations = 30;
  auto tmat = eigvals_lanczos(ops, block, n_eig_to_converge, 1e-12,
                              max_iterations, 1e-7, seed);
  toc();

  ofile["Alphas"] = tmat.alphas;
  ofile["Betas"] = tmat.betas;
  ofile["Eigenvalues"] = tmat.eigenvalues;
  ofile["Dimension"] = block.size();

  return EXIT_SUCCESS;
}
