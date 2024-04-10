#include <xdiag/all.hpp>

int main() {
  using namespace xdiag;
  using namespace arma;
  using fmt::format;

  say_hello();

  int N = 28;
  double J1 = 1.0 * 4;
  double Hx = 0.1 * 2;
  double J2 = -0.7 * 4;
  double H1 = (40.0 + 0.05) * 2;
  double H2 = (-40.0 + 0.05) * 2;

  double precision = 1e-3;
  double tau = 0.1;
  int n_steps = 1;

  cx_mat sx(mat({{0., 0.5}, {0.5, 0.}}), mat({{0., 0.}, {0., 0.}}));

  Log.set_verbosity(2);

  std::string filename;
  if ((N % 4) == 0) {
    filename = format("chain.{}.TFI.J1J2.HZHX.4sl.toml", N);
  } else if ((N % 4) == 2) {
    filename = format("chain.{}.TFI.J1J2.HZHX.2sl.toml", N);
  }
  auto lfile = FileToml(filename, 'r');

  // Create the two hamiltonians
  auto bonds1 = BondList(lfile["Interactions"]);
  bonds1["SX"] = sx;
  bonds1["J1"] = J1;
  bonds1["J2"] = J2;
  bonds1["HX"] = Hx;
  bonds1["HZ"] = H1;

  auto bonds2 = BondList(lfile["Interactions"]);
  bonds2["SX"] = sx;
  bonds2["J1"] = J1;
  bonds2["J2"] = J2;
  bonds2["HX"] = Hx;
  bonds2["HZ"] = H2;

  auto group = PermutationGroup(lfile["Symmetries"]);
  auto irrep = Representation(lfile["Gamma.C2.A"]);

  Block block;
  if ((N % 4) == 0) {
    block = Spinhalf(N, group, irrep, 4);
  } else if ((N % 4) == 2) {
    block = Spinhalf(N, group, irrep, 2);
  }
  XDiagPrint(block);

  // Create magnetization operator
  BondList mag;
  for (int i = 0; i < N; ++i) {
    mag << Bond("SZ", i);
  }

  // Create all-up starting state
  auto v = State(block);
  v.vector()(0) = 1.0;
  XDiagPrint(inner(mag, v));

  // Estimate the operator norms
  double anorm1 = norm_estimate(bonds1, block);
  double anorm2 = norm_estimate(bonds2, block);
  Log("Estimated operator norms:");
  Log("H1: {}", anorm1);
  Log("H2: {}", anorm2);

  for (int i = 0; i < n_steps; ++i) {

    Log("At time evolution step {}", i + 1);
    Log("Evolving with bonds1 {}", i + 1);
    v = time_evolve(bonds1, v, tau, precision, 10, anorm1);
    Log("Evolving with bonds2 {}", i + 1);
    v = time_evolve(bonds2, v, tau, precision, 10, anorm2);

    // Perform measurement
    XDiagPrint(inner(mag, v));
  }

  return EXIT_SUCCESS;
}
