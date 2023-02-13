//
// Created by Luke Staszewski on 30.01.23.
//
#include "../../catch.hpp"
#include <hydra/all.h>

using namespace hydra;
using namespace std;
using namespace arma;

TEST_CASE("analytic_case_free_particle_1D", "[time_evolution]") {
  // some constants
  int n_sites = 17;
  int nup = 1;
  int ndn = 0;
  double t = 1.5;
  const double pi = 3.14159265358979323846;
  const cx_double i = {0, 1};

  auto psi_analytic = [&n_sites, &t, &pi, &i](double time) {
    // returns time evolved 1D free particle with nearest neighbour hops

    cx_vec psi(n_sites);
    for (int m = 0; m < n_sites; m++) {
      cx_double c_m = 0;
      for (int k = 0; k < n_sites; k++) {
        // 1/N^2 * sum( e^(a_k))
        cx_double a_k = 2 * pi * i * k * m / n_sites;
        a_k += 2 * i * t * cos(2 * pi * k / n_sites) * time;
        c_m += exp(a_k);
      }
      c_m /= n_sites;
      psi(m) = c_m;
    }
    return psi;
  };

  // time evolve using hydra (numerically)
  auto block = Electron(n_sites, nup, ndn);
  BondList bonds;
  for (int i = 0; i < n_sites; i++) {
    bonds << Bond("HOP", "t", {i, (i + 1) % n_sites});
  }
  bonds["t"] = t;
  bonds["U"] = 1;

  auto psi_0_list = ProductState();
  for (int i = 0; i < n_sites; i++)
    psi_0_list << "Emp";
  psi_0_list[0] = "Up";
  auto psi_0 = State(block, psi_0_list);

  // seems to be accurate for times up towards 10^5
  arma::vec times = arma::logspace(-3, 2, 10);

  for (auto time : times) {
    double tol = 1e-2;
    auto w_expokit = time_evolve(bonds, psi_0, time, tol);
    arma::cx_vec w_analytic = psi_analytic(time);

    // norm is one so no division here by norm of true
    auto eps = arma::norm(w_expokit.vector() - w_analytic);

    cout << "err: " << eps / time << endl;
    cout << "time = " << time << endl;
    w_analytic.print("ana");
    w_expokit.vector().print("lanc");
    CHECK(eps / time < tol);
  }
}

TEST_CASE("analytic_case_free_particle_2D", "[time_evolution]") {
  // some constants
  int L = 3; // width of square lattice
  int n_sites = L * L;
  int nup = 1;
  int ndn = 0;
  double t = 1.5;
  double t1 = 0.3;
  const double pi = 3.14159265358979323846;
  const cx_double i = {0, 1};

  auto psi_analytic = [&L, &n_sites, &t, &t1, &pi, &i](double time) {
    // returns time evolved 1D free particle with nearest neighbour hops

    cx_vec psi(n_sites);
    for (int m = 0; m < L; m++) {
      for (int n = 0; n < L; n++) {
        cx_double c_mn = 0; // wave-function in position basis
        for (int k1 = 0; k1 < L; k1++) {
          for (int k2 = 0; k2 < L; k2++) {
            // nearest neighbour contribution
            double energy;
            energy = -2 * t * (cos(2 * pi * k1 / L) + cos(2 * pi * k2 / L));

            // next nearest neighbours
            energy +=
                -2 * t1 *
                (cos(2 * pi * (k1 + k2) / L) + cos(2 * pi * (k1 - k2) / L));

            // exponent
            cx_double a_k1k2 =
                -i * time * energy + 2 * pi * i * (k1 * m + k2 * n) / L;
            c_mn += exp(a_k1k2);
          }
        }
        int s = L * m + n;
        psi(s) = c_mn;
      }
    }

    psi /= pow(L, 2);
    return psi;
  };

  // time evolve using hydra (numerically)
  auto block = Electron(n_sites, nup, ndn);
  BondList bonds;
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      int s = L * i + j;
      // nearest neighbour hops
      int s_up_hop = L * ((i + 1) % L) + j;
      int s_dn_hop = L * i + (j + 1) % L;
      bonds << Bond("HOP", "t", {s, s_up_hop});
      bonds << Bond("HOP", "t", {s, s_dn_hop});

      // next nearest neighbour hops
      int s_across_up = L * ((i + 1) % L) + (j + 1) % L;
      int s_across_dn = L * ((L + (i - 1) % L) % L) + (j + 1) % L;
      bonds << Bond("HOP", "t1", {s, s_across_up});
      bonds << Bond("HOP", "t1", {s, s_across_dn});
    }
  }
  bonds["t"] = t;
  bonds["t1"] = t1;
  bonds["U"] = 1;

  auto psi_0_list = ProductState();
  for (int i = 0; i < n_sites; i++)
    psi_0_list << "Emp";
  psi_0_list[0] = "Up";
  auto psi_0 = State(block, psi_0_list);

  // seems to be accurate for times up towards 10^5
  arma::vec times = arma::logspace(-3, 2, 10);

  for (auto time : times) {
    double tol = 1e-10;
    auto w_expokit = time_evolve(bonds, psi_0, time, tol);
    arma::cx_vec w_analytic = psi_analytic(time);

    // norm is one so no division here by norm of true
    auto eps = arma::norm(w_expokit.vector() - w_analytic);

    cout << "err: " << eps / time << endl;
    cout << "time = " << time << endl;
    w_analytic.print("ana");
    w_expokit.vector().print("lanc");

    cout << "norm hydra " << norm(w_expokit) << endl;
    cout << "norm analytic " << norm(w_analytic) << endl;

    CHECK(eps / time < tol);
  }
}
