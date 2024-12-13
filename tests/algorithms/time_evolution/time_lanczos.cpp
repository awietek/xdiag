//
// Created by Luke Staszewski on 08.02.23.
//
#include "armadillo"
#include "lanczos_time_evolve.hpp"
#include <chrono>

using namespace std;
using namespace arma;
using namespace xdiag;

int main() {
  cout << "timing the lanczos algorithm ... ";

  // defining a model with xdiag
  int n_sites, nup, ndn;
  n_sites = 10;
  nup = 3;
  ndn = 3;
  double t, t1, U;
  t = 1;
  t1 = .5;
  U = 4;

  auto block = Electron(n_sites, nup, ndn);
  OpSum ops;
  for (int i = 0; i < n_sites; i++) {
    ops += "t" * Op("Hop", {i, (i + 1) % n_sites});
    ops += "t1" * Op("Hop", {i, (i + 2) % n_sites});
  }
  ops += "U" * Op("HubbardU");
  ops["t"] = t;
  ops["t1"] = t1;
  ops["U"] = U;

  auto psi_0_list = ProductState(
      {"Up", "Up", "Emp", "UpDn", "Emp", "Dn", "Dn", "Emp", "Emp", "Emp"});
  auto psi_0 = State(block, psi_0_list);

  // doing the time evolution
  double time = 30;
  auto tic = chrono::high_resolution_clock::now();
  auto psi = zahexpv(time, ops, psi_0, 1e-2, 10);
  auto toc = chrono::high_resolution_clock::now();
  cout << norm(psi.vector()) << endl;
  cout << "time taken: " << (toc - tic).count() * 1e-9 << endl;

  // timing with last xdiag

  auto tic1 = chrono::high_resolution_clock::now();
  auto psi_1 = time_evolve(ops, psi_0, time, 1e-2);
  auto toc1 = chrono::high_resolution_clock::now();
  cout << norm(psi.vector()) << endl;
  cout << "time taken xdiag: " << (toc1 - tic1).count() * 1e-9 << endl;
  return 0;
}
