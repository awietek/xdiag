#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {
  
  int n_sites = 16;
  int nup = n_sites / 2;
  Spinhalf block(n_sites, nup);
  XDIAG_PRINT(block);
  // Define the nearest-neighbor Heisenberg model
  BondList bonds;
  for (int i = 0; i < n_sites; ++i) {
    bonds += Bond("HB", "J", {i, (i + 1) % n_sites});
  }
  bonds["J"] = 1.0;

  set_verbosity(2);                  // set verbosity for monitoring progress
  double e0 = eigval0(bonds, block); // compute ground state energy
  
  Log("Ground state energy: {:.12f}", e0);
  
} catch (std::exception const &e) {
  traceback(e);
}
