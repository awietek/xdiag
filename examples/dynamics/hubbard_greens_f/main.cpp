#include <xdiag/all.hpp>

using namespace xdiag;
using namespace arma;
using fmt::format;
using namespace std::complex_literals;

int main() try {
  say_hello();

  // IO
  std::string latticeInput =
      XDIAG_DIRECTORY "/misc/data/square.8.hubbard.ttprime.toml";
  auto lfile = FileToml(latticeInput);

  std::string filename =
      XDIAG_DIRECTORY "/misc/data/examples_output/hubbard_greens_f.h5";
  auto outfile = FileH5(filename, "w!");

  // Define the Hubabrd model
  int N = 8;
  int nup = N / 2;
  int ndn = N / 2;
  auto ops = read_opsum(lfile, "Interactions");
  ops["Ty"] = 1.;
  ops["Tx"] = 1.;
  ops["Tprime+"] = 0.3;
  ops["Tprime-"] = 0.3;
  ops += "U" * Op("HubbardU");
  ops["U"] = -6.;
  auto opsTBC = ops;

  auto irrep = read_representation(lfile, "Gamma.C1.A");

  // compute groundstate (known to be at k=0)
  Log("Computing ground state ...");
  auto block = Electron(N, nup, ndn, irrep);
  auto [e0, gs] = eig0(ops, block);
  gs.make_complex();
  Log("done.");
  Log("Ground state energy: {:.12f}", e0);
  outfile["e0"] = e0;

  // loop through momenta
  std::vector<std::pair<std::string, arma::vec>> irreps = {
      {"Gamma.C1.A", {0., 0.}},
      {"M.C1.A", {3.1415926535897931, 3.1415926535897931}},
      {"Sigma0.C1.A", {1.5707963267948966, 1.5707963267948966}},
      {"Sigma1.C1.A", {1.5707963267948966, -1.5707963267948966}},
      {"Sigma2.C1.A", {-1.5707963267948966, 1.5707963267948966}},
      {"Sigma3.C1.A", {-1.5707963267948966, -1.5707963267948966}},
      {"X0.C1.A", {3.1415926535897931, 0.0000000000000000}},
      {"X1.C1.A", {0.0000000000000000, 3.1415926535897931}}};
  for (const auto &[name, momentum] : irreps) {
    Log("considering irrep {}, momentum [{}, {}]", name, momentum(0),
        momentum(1));
    auto aq_irrep = read_representation(lfile, name);
    auto c_q = symmetrize(Op("Cup", 0), aq_irrep);
    auto Av = apply(c_q, gs);
    auto nrm = norm(Av);
    Av /= nrm;
    XDIAG_SHOW(nrm);

    auto res = eigvals_lanczos_inplace(ops, Av);
    XDIAG_SHOW(res.eigenvalues(0));
    outfile[format("{}/norm", name)] = nrm;
    outfile[format("{}/alphas", name)] = res.alphas;
    outfile[format("{}/betas", name)] = res.betas;
    outfile[format("{}/eigs", name)] = res.eigenvalues;

    // compare to twisted boundary condition calculation,
    // cf. Zemljic & Prelovsek, PRB 75, 104514 (2007);
    // Tohyama, PRB 70, 174517 (2004).
    opsTBC["Tx"] = exp(1i * momentum(0));
    opsTBC["Ty"] = exp(1i * momentum(1));
    opsTBC["Tprime+"] =
        0.3 * exp(1i * (momentum(0) + momentum(1))); // r_ij = (1,1)
    opsTBC["Tprime-"] =
        0.3 * exp(1i * (momentum(0) - momentum(1))); // r_ij = (1,-1)
    // XDIAG_SHOW(opsTBC);
    // auto [e0_tbc, gs_tbc] = eig0(opsTBC, block);
    // XDIAG_SHOW(e0_tbc);
    auto c_q_tbc = symmetrize(Op("Cup", 0), irrep);
    Av = apply(c_q_tbc, gs);
    nrm = norm(Av);
    XDIAG_SHOW(nrm);
    Av /= nrm;

    auto res_tbc = eigvals_lanczos_inplace(opsTBC, Av);
    XDIAG_SHOW(res_tbc.eigenvalues(0));
  }

} catch (Error e) {
  error_trace(e);
}
