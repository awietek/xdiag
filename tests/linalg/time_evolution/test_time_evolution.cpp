// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

//
// Created by Luke Staszewski on 30.01.23.
//
#include <tests/catch.hpp>

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/kernels/sparse/csr_matrix.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/linalg/time_evolution/evolve_lanczos.hpp>
#include <xdiag/linalg/time_evolution/expm.hpp>
#include <xdiag/linalg/time_evolution/imaginary_time_evolve.hpp>
#include <xdiag/linalg/time_evolution/time_evolve.hpp>
#include <xdiag/linalg/time_evolution/time_evolve_expokit.hpp>
#include <xdiag/linalg/time_evolution/zahexpv.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/norm.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/timing.hpp>

using namespace xdiag;
using namespace std;
using namespace arma;

TEST_CASE("analytic_case_free_particle_1D", "[time_evolution]") try {
  Log("testing time evolution: analytic_case_free_particle_1D");
  Log.set_verbosity(0);

  // some constants
  int nsites = 17;
  int nup = 1;
  int ndn = 0;
  double t = 1.5;
  const double pi = 3.14159265358979323846;
  const cx_double i = {0, 1};

  auto psi_analytic = [&nsites, &t, &pi, &i](double time) {
    // returns time evolved 1D free particle with nearest neighbour hops

    cx_vec psi(nsites);
    for (int m = 0; m < nsites; m++) {
      cx_double c_m = 0;
      for (int k = 0; k < nsites; k++) {
        // 1/N^2 * sum( e^(a_k))
        cx_double a_k = (2 * k * m * pi * i) / (double)nsites;
        a_k += 2.0 * i * t * cos(2 * pi * k / nsites) * time;
        c_m += exp(a_k);
      }
      c_m /= nsites;
      psi(m) = c_m;
    }
    return psi;
  };

  // time evolve using xdiag (numerically)
  auto block = Electron(nsites, nup, ndn);
  OpSum ops;
  for (int i = 0; i < nsites; i++) {
    ops += "t" * Op("Hop", {i, (i + 1) % nsites});
  }
  ops["t"] = t;
  ops["U"] = 1;
  auto csr = csr_matrix(ops, block);

  std::vector<int> psi_0_list;
  for (int i = 0; i < nsites; i++) {
    psi_0_list.push_back(0);
  }
  psi_0_list[0] = 1;
  auto psi_0 = product_state(block, psi_0_list);

  // seems to be accurate for times up towards 10^5
  arma::vec times = arma::logspace(-3, 1, 10);

  for (auto time : times) {
    std::vector<double> tols = {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12};
    for (auto tol : tols) {
      auto w_expokit = time_evolve(ops, psi_0, time, tol, "expokit");
      auto w_expokit_csr = time_evolve(csr, psi_0, time, tol, "expokit");

      arma::cx_vec w_analytic = psi_analytic(time);

      // norm is one so no division here by norm of true
      auto eps = arma::norm(w_expokit.vectorC() - w_analytic);
      auto epscsr = arma::norm(w_expokit_csr.vectorC() - w_analytic);

      // cout << "err: " << eps / time << endl;
      // cout << "time = " << time << endl;
      // w_analytic.print("ana");
      // w_expokit.vector().print("lanc");

      REQUIRE(eps < 4 * tol);
      REQUIRE(epscsr < 4 * tol);
    }
  }
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
  throw;
}

TEST_CASE("analytic_case_free_particle_2D", "[time_evolution]") try {
  Log("testing time evolution: analytic_case_free_particle_2D");
  Log.set_verbosity(0);

  // some constants
  int L = 3; // width of square lattice
  int nsites = L * L;
  int nup = 1;
  int ndn = 0;
  double t = 1.5;
  double t1 = 0.3;
  const double pi = 3.14159265358979323846;
  const cx_double i = {0, 1};

  auto psi_analytic = [&L, &nsites, &t, &t1, &pi, &i](double time) {
    // returns time evolved 1D free particle with nearest neighbour hops

    cx_vec psi(nsites);
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
            cx_double a_k1k2 = -i * time * energy +
                               2 * pi * i * double(k1 * m + k2 * n) / (double)L;
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

  // time evolve using xdiag (numerically)
  auto block = Electron(nsites, nup, ndn);
  OpSum ops;
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      int s = L * i + j;
      // nearest neighbour hops
      int s_up_hop = L * ((i + 1) % L) + j;
      int s_dn_hop = L * i + (j + 1) % L;
      ops += "t" * Op("Hop", {s, s_up_hop});
      ops += "t" * Op("Hop", {s, s_dn_hop});

      // next nearest neighbour hops
      int s_across_up = L * ((i + 1) % L) + (j + 1) % L;
      int s_across_dn = L * ((L + (i - 1) % L) % L) + (j + 1) % L;
      ops += "t1" * Op("Hop", {s, s_across_up});
      ops += "t1" * Op("Hop", {s, s_across_dn});
    }
  }
  ops += "U" * Op("HubbardU");
  ops["t"] = t;
  ops["t1"] = t1;
  ops["U"] = 1;
  auto csr = csr_matrix(ops, block);

  std::vector<int> psi_0_list;
  for (int i = 0; i < nsites; i++) {
    psi_0_list.push_back(0);
  }
  psi_0_list[0] = 1;
  auto psi_0 = product_state(block, psi_0_list);

  // seems to be accurate for times up towards 10^5
  arma::vec times = arma::logspace(-3, 1, 10);

  for (auto time : times) {
    std::vector<double> tols = {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12};
    for (auto tol : tols) {
      auto w_expokit = time_evolve(ops, psi_0, time, tol, "expokit");
      auto w_expokit_csr = time_evolve(csr, psi_0, time, tol, "expokit");
      arma::cx_vec w_analytic = psi_analytic(time);

      // norm is one so no division here by norm of true
      auto eps = arma::norm(w_expokit.vectorC() - w_analytic);
      auto epscsr = arma::norm(w_expokit.vectorC() - w_analytic);

      // cout << "tol: " << tol << endl;
      // cout << "err: " << eps << endl;
      // cout << "time = " << time << endl;
      // w_analytic.print("ana");
      // w_expokit.vector().print("lanc");

      // cout << "norm xdiag " << norm(w_expokit) << endl;
      // cout << "norm analytic " << norm(w_analytic) << endl;

      REQUIRE(eps < 4 * tol);
      REQUIRE(epscsr < 4 * tol);
    }
  }
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
  throw;
}

TEST_CASE("tj_complex_timeevo", "[time_evolution]") try {
  Log("testing time evolution: tj_complex_timeevo");
  Log.set_verbosity(0);
  int L = 3;
  int nsites = L * L;

  // Create square lattice t-J model
  OpSum ops;
  for (int x = 0; x < L; ++x) {
    for (int y = 0; y < L; ++y) {
      int nx = (x + 1) % L;
      int ny = (y + 1) % L;

      int site = y * L + x;
      int right = y * L + nx;
      int top = ny * L + x;
      ops += "T" * Op("Hop", {site, right});
      ops += "Tasym" * Op("HopAsym", {site, right});
      ops += "J" * Op("tJSzSz", {site, right});
      ops += "T" * Op("Hop", {site, top});
      ops += "Tasym" * Op("HopAsym", {site, top});
      ops += "J" * Op("tJSzSz", {site, top});
    }
  }
  ops["T"] = 1.0;
  ops["Tasym"] = std::complex<double>(0.0, 0.2);
  ops["J"] = 0.4;

  // Create initial state
  auto pstate = ProductState();
  for (int x = 0; x < L; ++x) {
    for (int y = 0; y < L; ++y) {
      if (((x + y) % 2) == 0) {
        pstate.push_back(2);
      } else {
        pstate.push_back(1);
      }
    }
  }
  pstate[nsites / 2] = 0;
  auto block = tJ(nsites, nsites / 2, nsites / 2);
  auto csr = csr_matrixC(ops, block);
  auto H = matrixC(ops, block);
  double e0 = eigval0(ops, block);
  arma::cx_mat Hshift =
      H - e0 * arma::cx_mat(H.n_rows, H.n_cols, arma::fill::eye);

  auto psi_0 = State(block, false);
  xdiag::fill(psi_0, pstate);

  arma::vec times = arma::logspace(-1, 1, 3);

  // Expokit tests
  Log("testing time evolution: tj_complex_timeevo (EXPOKIT)");
  for (auto time : times) {
    std::vector<double> tols = {1e-2, 1e-6, 1e-10, 1e-12};
    for (auto tol : tols) {
      // Log("time: {}, tol: {}", time, tol);
      // tic();
      auto psi = time_evolve(ops, psi_0, time, tol, "expokit");
      // toc("OpSum");
      // tic();
      auto psicsr = time_evolve(csr, psi_0, time, tol, "expokit");
      // toc("CSRMatrix");

      cx_vec psi2 = expm(cx_mat(-1.0i * time * H)) * psi_0.vectorC();

      double eps = norm(psi2 - psi.vectorC());
      double epscsr = norm(psi2 - psicsr.vectorC());
      // cout << "tol: " << tol << endl;
      // cout << "err: " << eps << endl;
      // cout << "time = " << time << endl;
      // w_analytic.print("ana");
      // w_expokit.vector().print("lanc");

      // cout << "norm xdiag " << norm(w_expokit) << endl;
      // cout << "norm analytic " << norm(w_analytic) << endl;

      REQUIRE(eps < 4 * tol);
      REQUIRE(epscsr < 4 * tol);
    }
  }

  times = arma::logspace(-1, 0, 2);

  // lanczos tests
  Log("testing time evolution: tj_complex_timeevo (Lanczos)");
  for (auto time : times) {
    std::vector<double> tols = {1e-2, 1e-6, 1e-10, 1e-12};
    for (auto tol : tols) {
      // Log("time: {}, tol: {}", time, tol);
      // tic();
      auto psi = time_evolve(ops, psi_0, time, tol, "lanczos");
      // toc("OpSum");
      // tic();
      auto psicsr = time_evolve(csr, psi_0, time, tol, "lanczos");
      // toc("CSRMatrix");

      cx_vec psi2 = expm(cx_mat(-1.0i * time * H)) * psi_0.vectorC();

      auto ipsi = imaginary_time_evolve(ops, psi_0, time, tol, e0);
      auto ipsicsr = imaginary_time_evolve(csr, psi_0, time, tol, e0);

      cx_vec ipsi2 = expm(cx_mat(-time * Hshift)) * psi_0.vectorC();

      double eps = norm(psi2 - psi.vectorC());
      double epscsr = norm(psi2 - psicsr.vectorC());
      double ieps = norm(ipsi2 - ipsi.vectorC());
      double iepscsr = norm(ipsi2 - ipsicsr.vectorC());

      // Log("eps: {}, ieps: {}", eps, ieps);
      // Log("nmat: {}, nlcz: {}", arma::norm(ipsi.vectorC()),
      // arma::norm(ipsi2));

      // cout << "tol: " << tol << endl;
      // cout << "err: " << eps << endl;
      // cout << "time = " << time << endl;
      // w_analytic.print("ana");
      // w_expokit.vector().print("lanc");

      // cout << "norm xdiag " << norm(w_expokit) << endl;
      // cout << "norm analytic " << norm(w_analytic) << endl;

      REQUIRE(eps < 40 * tol);
      REQUIRE(ieps < 40 * tol);
      REQUIRE(epscsr < 40 * tol);
      REQUIRE(iepscsr < 40 * tol);
    }
  }

} catch (xdiag::Error e) {
  xdiag::error_trace(e);
  throw;
}

// TEST_CASE("zero_state_timeevo", "[time_evolution]") try {
//   int N = 4;
//   auto b = Spinhalf(N);
//   OpSum ops;
//   for (int i = 0; i < N; ++i) {
//     ops += Op("SdotS", {i, (i + 1) % N});
//   }
//   auto z = zero_state(b);
//   // auto zz1 = time_evolve(ops, z, 1.0);
//   // auto zz2 = time_evolve_expokit(ops, z, 1.0);
//   // auto zz3 = evolve_lanczos(ops, z, 1.0);
//   // auto zz4 = evolve_lanczos(ops, z, xdiag::complex(0.0, 1.0));

// } catch (xdiag::Error e) {
//   xdiag::error_trace(e);
// }

// A real Hamiltonian acting on a real state with a real time runs the purely
// real Lanczos exponential (the common imaginary-time case), which the complex
// tests above never reach. Every OpSum / CSR / in-place overload is checked
// against the dense matrix exponential.
TEST_CASE("time_evolution_real_arithmetic", "[time_evolution]") try {
  Log("testing time evolution: real Hamiltonian, real state, real time");
  int64_t nsites = 8;
  auto block = Spinhalf(nsites, 4);
  int64_t D = size(block);
  OpSum ops;
  for (int64_t i = 0; i < nsites; ++i) {
    ops += 1.0 * Op("SdotS", {i, (i + 1) % nsites});
    ops += 0.4 * Op("SdotS", {i, (i + 2) % nsites});
  }
  auto csr = csr_matrix(ops, block);
  mat H = matrix(ops, block);
  State psi_0 = random_state(block, true, 1, 7);
  REQUIRE(isreal(psi_0));
  vec v0 = psi_0.vector();
  double shift = -2.5;
  double tol = 1e-10;

  for (double time : {0.1, 1.0, 3.0}) {
    // imaginary time evolution exp(-time (H - shift))
    vec ref_imag = expmat(mat(-time * (H - shift * eye(D, D)))) * v0;
    State psi = imaginary_time_evolve(ops, psi_0, time, tol, shift);
    REQUIRE(isreal(psi));
    REQUIRE(norm(psi.vector() - ref_imag) < 1e-7 * norm(ref_imag));
    State psicsr = imaginary_time_evolve(csr, psi_0, time, tol, shift);
    REQUIRE(isreal(psicsr));
    REQUIRE(isapprox(psicsr, psi, 1e-7, 1e-7));

    State psi_in = psi_0;
    imaginary_time_evolve_inplace(ops, psi_in, time, tol, shift);
    REQUIRE(isreal(psi_in));
    REQUIRE(isapprox(psi_in, psi, 1e-7, 1e-7));
    State psi_in_csr = psi_0;
    imaginary_time_evolve_inplace(csr, psi_in_csr, time, tol, shift);
    REQUIRE(isapprox(psi_in_csr, psi, 1e-7, 1e-7));

    // evolve_lanczos with a real tau: OpSum, CSR, in-place
    vec ref_tau = expmat(mat(-time * H)) * v0;
    EvolveLanczosResult r = evolve_lanczos(ops, psi_0, -time, tol);
    REQUIRE(isreal(r.state));
    REQUIRE(r.niterations > 0);
    REQUIRE(norm(r.state.vector() - ref_tau) < 1e-7 * norm(ref_tau));
    EvolveLanczosResult rcsr = evolve_lanczos(csr, psi_0, -time, tol);
    REQUIRE(isapprox(rcsr.state, r.state, 1e-7, 1e-7));
    State s_in = psi_0;
    EvolveLanczosInplaceResult ri = evolve_lanczos_inplace(csr, s_in, -time, tol);
    REQUIRE(isapprox(s_in, r.state, 1e-7, 1e-7));
    REQUIRE(ri.niterations == rcsr.niterations);

    // a complex tau on a real state promotes the state to complex
    cx_vec v0c = conv_to<cx_vec>::from(v0);
    cx_vec ref_cx = expmat(cx_mat(xdiag::complex(0.0, -time) * H)) * v0c;
    EvolveLanczosResult rc = evolve_lanczos(ops, psi_0, xdiag::complex(0.0, -time), tol);
    REQUIRE(!isreal(rc.state));
    REQUIRE(norm(rc.state.vectorC() - ref_cx) < 1e-7);
    EvolveLanczosResult rccsr =
        evolve_lanczos(csr, psi_0, xdiag::complex(0.0, -time), tol);
    REQUIRE(isapprox(rccsr.state, rc.state, 1e-7, 1e-7));
    State c_in = psi_0;
    evolve_lanczos_inplace(ops, c_in, xdiag::complex(0.0, -time), tol);
    REQUIRE(isapprox(c_in, rc.state, 1e-7, 1e-7));

    // normalize = true returns a unit vector
    EvolveLanczosResult rn = evolve_lanczos(ops, psi_0, -time, tol, 0.0, true);
    REQUIRE(std::abs(norm(rn.state) - 1.0) < 1e-10);
  }

  // argument checks: zero state, non-hermitian operator
  {
    State zero = zero_state(block);
    REQUIRE_THROWS_AS(evolve_lanczos(ops, zero, 1.0), Error);
    REQUIRE_THROWS_AS(evolve_lanczos(ops, zero, xdiag::complex(0.0, 1.0)),
                      Error);

    OpSum nonherm = 1.0 * (Op("S+", 0) * Op("S-", 1));
    REQUIRE_THROWS_AS(evolve_lanczos(nonherm, psi_0, 1.0), Error);
    REQUIRE_THROWS_AS(
        evolve_lanczos(nonherm, psi_0, xdiag::complex(0.0, 1.0)), Error);
    REQUIRE_THROWS_AS(imaginary_time_evolve(nonherm, psi_0, 1.0), Error);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// A small Krylov space forces the Expokit step-size controller to shorten and
// reject steps; the result must still reach the requested accuracy.
TEST_CASE("zahexpv_small_krylov", "[time_evolution]") try {
  Log("testing zahexpv with a small Krylov space");
  arma::arma_rng::set_seed(5);
  int64_t n = 60;
  cx_mat B(n, n, fill::randn);
  cx_mat A = B + B.t(); // hermitian
  double anorm = norm(A, "inf");
  auto apply_A = [&A](cx_vec const &v) {
    return cx_vec(xdiag::complex(0.0, -1.0) * (A * v));
  };
  auto dot_f = [](cx_vec const &v, cx_vec const &w) { return cdot(v, w); };
  cx_vec w0(n, fill::randn);
  w0 /= norm(w0);

  // reference exp(-i t A) w0 from the eigendecomposition of A
  vec evals;
  cx_mat evecs;
  eig_sym(evals, evecs, A);
  cx_vec w0_eig = evecs.t() * w0;

  for (int m : {6, 8, 30}) {
    for (double time : {0.5, 5.0}) {
      cx_vec w = w0;
      double tol = 1e-10;
      auto [err, hump] = zahexpv(time, apply_A, dot_f, w, anorm, tol / time, m);
      cx_vec ref = evecs * (exp(xdiag::complex(0.0, -time) * evals) % w0_eig);
      REQUIRE(norm(w - ref) < 1e-6);
      REQUIRE(std::abs(norm(w) - 1.0) < 1e-8);
    }
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
