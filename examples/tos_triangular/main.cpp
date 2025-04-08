#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {

  std::vector<double> energies;
  std::vector<std::string> IrrepList;
  std::vector<int> Szsector;
  int Nsites = 18;
  int numbeig = 6; // number of Eigenvalues to converge.

  std::vector<std::string> Irreps = {
      "Gamma.C2.A", "Gamma.C2.B", "K0.C1.A", "K1.C1.A", "M.C2.A",
      "M.C2.B",     "X0.C1.A",    "X1.C1.A", "X2.C1.A", "Z0.C1.A",
      "Z1.C1.A",    "0.C1.A",     "1.C1.A",  "2.C1.A"};

  auto fl = FileToml("triangular.18.10418.J1J2.sublattices.tsl.toml");

  OpSum ops = read_opsum(fl, "Interactions");

  // For this ratio, we are in the 120 phase
  ops["J1"] = 1.0;
  ops["J2"] = 0.0;

  for (auto irrep : Irreps) {

    auto irrep2 = read_representation(fl, irrep, "Symmetries");

    for (int nup = 0; nup <= Nsites; nup++) {
      auto block = Spinhalf(Nsites, nup, irrep2);

      auto res = eigvals_lanczos(ops, block, numbeig);
      arma::vec eig0 = res.eigenvalues;

      for (int i = 0; i < eig0.n_elem; i++) {
        energies.push_back(eig0[i]);
        IrrepList.push_back(irrep);
        Szsector.push_back(nup);
      }
    }
  }

  // Construct the filename
  std::string flstring = "energies_tower_of_states.triangular.Nsites." +
                         std::to_string(Nsites) + ".outfile.txt";
  std::ofstream outfile(flstring);
  for (int i = 0; i < energies.size(); i++) {
    outfile << energies[i] << "," << Szsector[i] << "," << IrrepList[i] << "\n";
  }
  outfile.close();
  return 0;
} catch (Error e) {
  error_trace(e);
}
