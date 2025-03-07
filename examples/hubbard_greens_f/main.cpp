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
  int nup = 3;
  int ndn = 3;
  auto ops = read_opsum(lfile, "Interactions");
  auto t = -1.;
  auto tp = 0.;
  ops["Ty"] = t;
  ops["Tx"] = t;
  ops["Tprime+"] = tp;
  ops["Tprime-"] = tp;
  ops += "U" * Op("HubbardU");
  ops["U"] = -6.;
  auto opsTBC = ops;

  auto irrep = read_representation(lfile, "Gamma.C1.A");

  // compute groundstate (known to be at k=0)
  Log("Computing ground state ...");
  auto block = Electron(N, nup, ndn, irrep);
  auto [e0, gs] = eig0(ops, block);
  gs.make_complex();
  Log("gag.");
  Log("Ground state energy: {:.12f}", e0);
  outfile["e0"] = e0;

  auto tbc_mesh =
      [&opsTBC, &t, &tp, &block, &N,
       &irrep](double const &kx,
               double const &ky) -> std::pair<double, EigvalsLanczosResult> {
    // lambda to perform TBC calc.
    // compare to twisted boundary condition calculation,
    // cf. Zemljic & Prelovsek, PRB 75, 104514 (2007);
    // Tohyama, PRB 70, 174517 (2004).
    Log("TBC calculation, momentum [{}, {}]", kx, ky);
    opsTBC["Tx"] = t * exp(1i * kx);
    opsTBC["Ty"] = t * exp(1i * ky);
    opsTBC["Tprime+"] = tp * exp(1i * (kx + ky)); // r_ij = (1,1)
    opsTBC["Tprime-"] = tp * exp(1i * (kx - ky)); // r_ij = (1,-1)
    XDIAG_SHOW(opsTBC);
    auto [e0_tbc, gs_tbc] = eig0(opsTBC, block, 1e-12, 100);
    // XDIAG_SHOW(e0_tbc);
    auto c_q_tbc = symmetrize(sqrt((double)N) * Op("Cup", 0), irrep);
    auto Av = apply(c_q_tbc, gs_tbc);
    auto nrm = norm(Av);
    XDIAG_SHOW(nrm);
    Av /= nrm;

    auto res_tbc = eigvals_lanczos_inplace(opsTBC, Av, 1, 1e-12, 100);
    XDIAG_SHOW(res_tbc.eigenvalues(0));
    Log("n iterations {}", res_tbc.niterations);
    return {nrm, res_tbc};
  };

  // loop through momenta of the cluster reciprocal lattice
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
    auto c_q = symmetrize(sqrt((double)N) * Op("Cup", 0),
                          aq_irrep); // normalize to 1/sqrt(N)
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

    // compare to tbc calc
    auto [nrm_tbc, res_tbc] = tbc_mesh(momentum(0), momentum(1));
    XDIAG_SHOW(nrm_tbc);
    XDIAG_SHOW(res_tbc.eigenvalues(0));
  }

  // TBC mesh, 10x10 grid
  int nkx = 101;
  int nky = nkx;
  auto momentum = [&nkx, &nky](int const &ikx,
                               int const &iky) -> std::pair<double, double> {
    return {(double)ikx * 3.1415926535897931 / (double)(nkx-1),
            (double)iky * 3.1415926535897931 / (double)(nky-1)};
  };
  for (int ikx = 0; ikx < nkx; ++ikx) {
    auto iky=ikx;
    // for (int iky = 0; iky < nky; ++iky) {
      auto [kx, ky] = momentum(ikx, iky);
      auto [nrm, res] = tbc_mesh(kx, ky);
      XDIAG_SHOW(nrm);
      XDIAG_SHOW(res.eigenvalues(0));

      outfile[format("{}/{}/norm", ikx, iky)] = nrm;
      outfile[format("{}/{}/alphas", ikx, iky)] = res.alphas;
      outfile[format("{}/{}/betas", ikx, iky)] = res.betas;
      outfile[format("{}/{}/eigs", ikx, iky)] = res.eigenvalues;
      outfile[format("{}/{}/kx", ikx, iky)] = kx;
      outfile[format("{}/{}/ky", ikx, iky)] = ky;
    // }
  }

} catch (Error e) {
  error_trace(e);
}