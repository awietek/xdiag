#include <xdiag/all.hpp>

using namespace xdiag;
using namespace arma;
using fmt::format;
using namespace std::complex_literals;

int main() try {
  say_hello();
  set_verbosity(1)
;
  // IO
  std::string latticeInput =
      XDIAG_DIRECTORY "/misc/data/square.20.tJ.toml";
  auto lfile = FileToml(latticeInput);

  std::string filename =
      XDIAG_DIRECTORY "/misc/data/examples_output/tJ_conductivity.h5";
  auto outfile = FileH5(filename, "w!");

  // Lanczos parameters
  auto precision = 0; // turns off convergence checking
  auto maxiters = 200;

  // Define the model
  int N = 20;
  int nup = 5;
  int ndn = 4; // one hole
  auto ops = read_opsum(lfile, "Interactions");
  auto t = 1.;
  auto J = 0.3;
  ops["T"] = t;
  ops["J"] = J;
  XDIAG_SHOW(ops);

  auto irrep = read_representation(lfile, "Gamma.C1.A");

  // compute groundstate
  Log("computing ground state");
  Log("Computing ground state ...");
  auto block = tJ(N, nup, ndn, irrep);
  XDIAG_SHOW(block);
  auto [e0, gs] = eig0(ops, block);
  Log("Ground state energy: {:.12f}", e0);
  outfile["e0"] = e0;

  // create & apply current operator
  Log("applying curr. operator");
  auto current = symmetrize(1i * Op("Hop", {0, 1}),
                            irrep);
  // XDIAG_SHOW(current);
  auto Av = apply(current, gs);
  auto nrm = norm(Av);
  Av /= nrm;
  XDIAG_SHOW(nrm);
  outfile["norm"] = nrm;

  // second run
  Log("second Laczos run on A|gs>");
  auto res = eigvals_lanczos_inplace(ops, Av, 1, precision, maxiters);
  XDIAG_SHOW(res.eigenvalues(0));
  outfile["alphas"] = res.alphas;
  outfile["betas"] = res.betas;
  outfile["eigs"] = res.eigenvalues;

} catch (Error e) {
  error_trace(e);
}