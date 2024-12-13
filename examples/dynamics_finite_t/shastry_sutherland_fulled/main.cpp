#include <xdiag/all.hpp>

int main(int argc, char **argv) {
  using namespace xdiag;
  using namespace arma;
  using fmt::format;

  say_hello();
  
  // Parse input arguments
  assert(argc == 7);
  int n_sites = atoi(argv[1]);              // number of sites
  int n_up = atoi(argv[2]);                 // number of upspins
  std::string kname = std::string(argv[3]); // momentum k
  std::string qname = std::string(argv[4]); // momentum q
  double J = atof(argv[5]);
  double Jd = atof(argv[6]);

  auto lfile = FileToml(format("shastry.{}.HB.J.Jd.fsl.toml", n_sites), 'r');
  auto ofile = FileH5(
      format(
          "outfiles/outfile.shastry.{}.J.{:.2f}.Jd.{:.2f}.nup.{}.k.{}.q.{}.h5",
          n_sites, J, Jd, n_up, kname, qname),
      "w!");

  auto ops = OpSum(lfile["Interactions"]);
  ops["J"] = J;
  ops["Jd"] = Jd;
  auto group = PermutationGroup(lfile["Symmetries"]);
  auto irrep_k = Representation(lfile[kname]);
  auto block_k = Spinhalf(n_sites, n_up, group, irrep_k);

  Log("Diagonalizing H in block nup: {}, k: {}", n_up, kname);
  auto H = matrix(ops, block_k);
  vec eigval;
  cx_mat eigvec;
  eig_sym(eigval, eigvec, H);
  ofile["EigenvaluesK"] = eigval;

  Log("Diagonalizing H in block nup: {}, k+q: {}+{}", n_up, kname, qname);
  auto irrep_q = Representation(lfile[qname]);
  auto irrep_k_q = irrep_k * irrep_q;
  auto block_k_q = Spinhalf(n_sites, n_up, group, irrep_k_q);

  auto H_k_q = matrix(ops, block_k_q);
  vec eigval_k_q;
  cx_mat eigvec_k_q;
  eig_sym(eigval_k_q, eigvec_k_q, H_k_q);
  ofile["EigenvaluesKQ"] = eigval_k_q;

  Log("Computing S(q), q: {}", qname);
  mat coords = lfile["Coordinates"].as<mat>();
  auto q = lfile[qname + std::string(".momentum")].as<vec>();
  OpSum S_of_q_ops;
  for (int site = 0; site < n_sites; ++site) {
    complex phase = exp(1i * (q(0) * coords(site, 0) + q(1) * coords(site, 1)));
    S_of_q_ops << Op("Sz", phase / n_sites, site);
  }
  auto S_of_q = matrix(S_of_q_ops, block_k, block_k_q);
  cx_mat S_of_q_eig = eigvec_k_q.t() * S_of_q * eigvec;
  ofile["SofQ"] = S_of_q_eig;

  return EXIT_SUCCESS;
}
