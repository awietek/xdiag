#include <hydra/all.h>

int main(int argc, char **argv) {
  using namespace hydra;
  using namespace arma;
  using fmt::format;

  // Parse input arguments
  assert(argc == 5);
  int n_sites = atoi(argv[1]);              // number of sites
  int n_up = atoi(argv[2]);                 // number of upspins
  std::string kname = std::string(argv[3]); // momentum k
  std::string qname = std::string(argv[4]); // momentum q

  auto lfile = FileToml(format("shastry.{}.HB.J.Jd.fsl.toml", n_sites), 'r');
  auto ofile = FileH5(
      format("outfiles/outfile.shastry.{}.HB.J.Jd.fsl.nup.{}.k.{}.q.{}.h5",
             n_sites, n_up, kname, qname));

  auto bonds = BondList(lfile["Interactions"]);
  auto group = PermutationGroup(lfile["Symmetries"]);
  auto irrep_k = Representation(lfile[kname]);
  auto block_k = Spinhalf(n_sites, n_up, group, irrep);

  Log("Diagonalizing H in block nup: {}, k: {}", n_up, kname);
  auto hamiltonian = matrix(bonds, block);
  vec eigval;
  cx_mat eigvec;
  eig_sym(eigval, eigvec, H);
  ofile["EigenvaluesK"] = eigval;

  Log("Diagonalizing H in block nup: {}, k+q: {}+{}", n_up, kname, qname);
  auto irrep_q = Representation(lfile[qname]);
  auto irrep_k_q = irrep_k * irrep_q;
  auto block_k_q = Spinhalf(n_sites, n_up, group, irrep);

  auto H_k_q = matrix(bonds, block_k_q);
  vec eigval_k_q;
  cx_mat eigvec_k_q;
  eig_sym(eigval_k_q, eigvec_k_q, H_k_q);
  ofile["EigenvaluesKQ"] = eigval_k_q;

  Log("Computing S(q), q: {}", qname);
  auto coords = lfile["Coordinates"].as<mat>();
  auto q = lfile[qname + std::string(".momentum")].as<vec>();
  BondList S_of_q_bonds;
  for (int site = 0; site < n_sites; ++site) {
    vec r = coords.row(site);
    complex phase = exp(1i * dot(q, r));
    S_of_q_bonds << Bond("SZ", phase / n_sites, site);
  }
  auto S_of_q = matrix(S_of_q_bonds, block_k, block_k_q);
  auto S_of_q_eig = eigvec_k_q.t() * S_of_q * eigvec;
  ofile["SofQ"] = S_of_q_eig;

  return EXIT_SUCCESS;
}
