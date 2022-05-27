#include <hydra/all.h>

int main() {
  using namespace hydra;

  // Define the all-to-all hopping and Heisenberg terms on 4 sites
  BondList bonds;
  bonds << Bond("HOP", "T01", {0, 1});
  bonds << Bond("HOP", "T02", {0, 2});
  bonds << Bond("HOP", "T03", {0, 3});
  bonds << Bond("HOP", "T12", {1, 2});
  bonds << Bond("HOP", "T13", {1, 3});
  bonds << Bond("HOP", "T23", {2, 3});
  bonds << Bond("HB", "J01", {0, 1});
  bonds << Bond("HB", "J02", {0, 2});
  bonds << Bond("HB", "J03", {0, 3});
  bonds << Bond("HB", "J12", {1, 2});
  bonds << Bond("HB", "J13", {1, 3});
  bonds << Bond("HB", "J23", {2, 3});

  // Set all the couplings to some value
  Couplings couplings;
  couplings["T01"] = complex(-6.6934277670e-1, -5.3885847860e-2);
  couplings["T02"] = complex(-2.9743014890e-2, -8.0010844280e-2);
  couplings["T03"] = complex(-8.9637147750e-1, 1.5197637250e+0);
  couplings["T12"] = complex(1.9142414960e-1, 7.7319093400e-1);
  couplings["T13"] = complex(-3.4208390920e-1, 7.6560791580e-1);
  couplings["T23"] = complex(-4.5034805430e-1, 2.1385531630e-1);
  couplings["J01"] = -6.9753683290e-1;
  couplings["J02"] = 1.0168476600e+0;
  couplings["J03"] = -8.8417711530e-1;
  couplings["J12"] = 3.5280741620e+0;
  couplings["J13"] = 1.2226585440e-1;
  couplings["J23"] = -2.5944271800e+0;

  // Define the number of lattice sites, number of upspins
  int n_sites = bonds.n_sites();
  int n_up = 2;
  int n_dn = 2;

  // Create the corresponding block of the Spinhalf Hilbertspace
  auto block = tJ(n_sites, n_up, n_dn);

  // Compute the full Hamiltonian matrix (complex, double precision)
  auto H = MatrixCplx(bonds, couplings, block, block);

  // Perform a full diagonalization
  auto eigs = lila::EigenvaluesSym(H);

  // Output the eigenvalues
  LilaPrint(eigs);

  /////
  // Alternatively, compute spectrum without number conservation

  // Create spinhalf block without particle numbers ...
  auto block_no_num = Spinhalf(n_sites);

  // ... create the Hamiltonian and compute eigenvalues
  auto H_no_num = MatrixCplx(bonds, couplings, block_no_num, block_no_num);
  auto eigs_no_num = lila::EigenvaluesSym(H_no_num);
  LilaPrint(eigs_no_num);

  return EXIT_SUCCESS;
}
