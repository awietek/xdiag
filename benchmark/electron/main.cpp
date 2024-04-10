#include <xdiag/all.hpp>

int main() {
  using namespace xdiag;
  using namespace arma;
  using fmt::format;

  Log.set_verbosity(2);
  
  int n_sites = 14;
  int nup = 6;
  int ndn = 6;
      
  BondList bonds;
  for (int s1 = 0; s1 < n_sites; ++s1) {
    int s2 = (s1+1) % n_sites;
    bonds << Bond("HOP", "T", {s1, s2});
    bonds << Bond("ISING", "JZ", {s1, s2});
    bonds << Bond("EXCHANGE", "JEX", {s1, s2});
  }
  bonds["T"] = 0.0;
  
  bonds["JZ"] = 1.0;
  bonds["JEX"] = 0.0;

  auto block = Electron(n_sites, nup, ndn);
  XDiagPrint(block);
  auto e0 = eig0(bonds, block);
  XDiagPrint(e0);

  return EXIT_SUCCESS;
}
