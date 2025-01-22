#include "../../catch.hpp"

#include <iostream>

#include "../spinhalf/testcases_spinhalf.hpp"
#include "../tj/testcases_tj.hpp"
#include "testcases_electron.hpp"

#include <xdiag/algebra/matrix.hpp>
#include <xdiag/utils/close.hpp>

using namespace xdiag;

void test_electron_np_no_np_matrix(int nsites, OpSum ops) try {

  auto block_full = Electron(nsites);
  auto H_full = matrixC(ops, block_full, block_full);
  REQUIRE(H_full.is_hermitian(1e-12));
  arma::Col<double> eigs_full;
  arma::eig_sym(eigs_full, H_full);

  std::vector<double> all_eigs;
  for (int nup = 0; nup <= nsites; ++nup)
    for (int ndn = 0; ndn <= nsites; ++ndn) {

      auto block = Electron(nsites, nup, ndn);
      auto H = matrixC(ops, block, block);
      REQUIRE(H.is_hermitian(1e-12));

      arma::Col<double> eigs;
      arma::eig_sym(eigs, H);

      for (auto eig : eigs)
        all_eigs.push_back(eig);
    }
  std::sort(all_eigs.begin(), all_eigs.end());
  REQUIRE(close(arma::Col(all_eigs), eigs_full));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

TEST_CASE("electron_matrix", "[electron]") try {
  using namespace xdiag::testcases::electron;

  OpSum ops;

  // Compare with matrix from Weisse and Fehske
  Log("electron_matrix: Hubbard Weisse & Fehske");
  int nsites = 4;
  int nup = 3;
  int ndn = 2;
  double t = 1.0;
  double U = 5.0;
  auto block = Electron(nsites, nup, ndn);

  for (int i = 0; i < nsites; ++i) {
    ops += "T" * Op("Hop", {i, (i + 1) % nsites});
  }
  ops += "U" * Op("HubbardU");
  ops["T"] = 1.0;
  ops["U"] = 5.0;
  auto H1 = matrix(ops, block, block);
  double tp = t;
  double tm = -t;
  double UU = U;
  double U2 = 2 * U;
  arma::Mat<double> H1_correct = {
      {U2, tm, 0., 0., tp, 0., tm, 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., tm, 0., 0., 0., 0., 0.}, // 1
      {tm, U2, tm, tm, 0., tp, 0., tm, 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., tm, 0., 0., 0., 0.}, // 2
      {0., tm, U2, 0., tm, 0., 0., 0., tm, 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., tm, 0., 0., 0.}, // 3
      {0., tm, 0., UU, tm, 0., 0., 0., 0., tm, 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., tm, 0., 0.}, // 4
      {tp, 0., tm, tm, UU, tm, 0., 0., 0., 0., tm, 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., tm, 0.}, // 5
      {0., tp, 0., 0., tm, UU, 0., 0., 0., 0., 0., tm,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., tm}, // 6
      {tm, 0., 0., 0., 0., 0., U2, tm, 0., 0., tp, 0.,
       tm, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, // 7
      {0., tm, 0., 0., 0., 0., tm, UU, tm, tm, 0., tp,
       0., tm, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, // 8
      {0., 0., tm, 0., 0., 0., 0., tm, UU, 0., tm, 0.,
       0., 0., tm, 0., 0., 0., 0., 0., 0., 0., 0., 0.}, // 9
      {0., 0., 0., tm, 0., 0., 0., tm, 0., U2, tm, 0.,
       0., 0., 0., tm, 0., 0., 0., 0., 0., 0., 0., 0.}, // 10
      {0., 0., 0., 0., tm, 0., tp, 0., tm, tm, U2, tm,
       0., 0., 0., 0., tm, 0., 0., 0., 0., 0., 0., 0.}, // 11
      {0., 0., 0., 0., 0., tm, 0., tp, 0., 0., tm, UU,
       0., 0., 0., 0., 0., tm, 0., 0., 0., 0., 0., 0.}, // 12

      {0., 0., 0., 0., 0., 0., tm, 0., 0., 0., 0., 0.,
       UU, tm, 0., 0., tp, 0., tm, 0., 0., 0., 0., 0.}, // 13
      {0., 0., 0., 0., 0., 0., 0., tm, 0., 0., 0., 0.,
       tm, U2, tm, tm, 0., tp, 0., tm, 0., 0., 0., 0.}, // 14
      {0., 0., 0., 0., 0., 0., 0., 0., tm, 0., 0., 0.,
       0., tm, UU, 0., tm, 0., 0., 0., tm, 0., 0., 0.}, // 15
      {0., 0., 0., 0., 0., 0., 0., 0., 0., tm, 0., 0.,
       0., tm, 0., U2, tm, 0., 0., 0., 0., tm, 0., 0.}, // 16
      {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., tm, 0.,
       tp, 0., tm, tm, UU, tm, 0., 0., 0., 0., tm, 0.}, // 17
      {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., tm,
       0., tp, 0., 0., tm, U2, 0., 0., 0., 0., 0., tm}, // 18
      {tm, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       tm, 0., 0., 0., 0., 0., UU, tm, 0., 0., tp, 0.}, // 19
      {0., tm, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., tm, 0., 0., 0., 0., tm, UU, tm, tm, 0., tp}, // 20
      {0., 0., tm, 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., tm, 0., 0., 0., 0., tm, U2, 0., tm, 0.}, // 21
      {0., 0., 0., tm, 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., tm, 0., 0., 0., tm, 0., UU, tm, 0.}, // 22
      {0., 0., 0., 0., tm, 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., tm, 0., tp, 0., tm, tm, U2, tm}, // 23
      {0., 0., 0., 0., 0., tm, 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., tm, 0., tp, 0., 0., tm, U2} // 24
  };

  REQUIRE(close(H1, H1_correct));
  for (int i = 0; i < 24; ++i)
    for (int j = 0; j < 24; ++j) {
      REQUIRE(close(std::real(H1(i, j)), H1_correct(i, j)));
    }

  //////////////////////////////////
  // Test two site exact solution
  {
    OpSum ops;
    ops += "T" * Op("Hop", {0, 1});
    ops += "U" * Op("HubbardU");
    auto block2 = Electron(2, 1, 1);
    for (int i = 0; i < 20; ++i) {
      double U = 1.234 * i;
      Log("electron_matrix: two-site exact solution test, U={}", U);
      ops["T"] = 1.0;
      ops["U"] = U;
      double e0_exact = 0.5 * (U - sqrt(U * U + 16));
      auto H = matrix(ops, block2, block2);
      REQUIRE(H.is_hermitian(1e-8));
      arma::Col<double> eigs;
      arma::eig_sym(eigs, H);
      double e0 = eigs(0);
      // printf("e0: %f, e0_exact: %f\n", e0, e0_exact);
      REQUIRE(close(e0_exact, e0));
    }
  }

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (real)
  for (int nsites = 3; nsites < 7; ++nsites) {
    Log("electron_matrix: free fermion random all-to-all test, (real), N={}",
        nsites);

    ops = freefermion_alltoall(nsites);

    // Create single particle matrix
    arma::Mat<double> Hs(nsites, nsites, arma::fill::zeros);
    for (auto [cpl, op] : ops) {
      REQUIRE(op.sites().size() == 2);
      int s1 = op[0];
      int s2 = op[1];
      auto name = cpl.string();
      Hs(s1, s2) = -ops[name].as<double>();
      Hs(s2, s1) = -ops[name].as<double>();
    }

    arma::vec seigs;
    arma::eig_sym(seigs, Hs);
    for (int nup = 0; nup <= nsites; ++nup)
      for (int ndn = 0; ndn <= nsites; ++ndn) {

        // Compute exact gs energy
        double e0_exact = 0;
        for (int i = 0; i < nup; ++i)
          e0_exact += seigs(i);
        for (int i = 0; i < ndn; ++i)
          e0_exact += seigs(i);

        auto block3 = Electron(nsites, nup, ndn);
        auto Hr = matrix(ops, block3, block3);
        REQUIRE(Hr.is_hermitian(1e-8));
        arma::vec eigsr;
        arma::eig_sym(eigsr, Hr);
        auto Hc = matrixC(ops, block3, block3);
        REQUIRE(Hc.is_hermitian(1e-8));
        arma::vec eigsc;
        arma::eig_sym(eigsc, Hc);
        REQUIRE(close(eigsr, eigsc));
        double e0 = eigsr(0);
        // printf("nup: %d, ndn: %d, dim: %d, e0: %f, e0_exact: %f\n", nup, ndn,
        //        (int)evecs.size(), e0, e0_exact);
        REQUIRE(close(e0_exact, e0));
      }
  }

  //////////////////////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (cplx, up/dn different)
  for (int nsites = 3; nsites < 7; ++nsites) {
    Log("electron_matrix: free fermion random all-to-all test, (cplx), N={}",
        nsites);
    ops = freefermion_alltoall_complex_updn(nsites);

    // Create single particle matrix for upspins
    arma::cx_mat Hs_up(nsites, nsites, arma::fill::zeros);
    for (auto [cpl, op] : ops) {
      if (op.type() == "Hopup") {
        assert(op.size() == 2);
        int s1 = op[0];
        int s2 = op[1];
        auto name = cpl.string();
        Hs_up(s1, s2) = -ops[name].as<complex>();
        Hs_up(s2, s1) = -conj(ops[name].as<complex>());
      }
    }
    arma::vec seigs_up;
    arma::eig_sym(seigs_up, Hs_up);

    // Create single particle matrix for dnspins
    arma::cx_mat Hs_dn(nsites, nsites, arma::fill::zeros);
    for (auto [cpl, op] : ops) {
      if (op.type() == "Hopdn") {
        assert(op.size() == 2);
        int s1 = op[0];
        int s2 = op[1];
        auto name = cpl.string();
        Hs_dn(s1, s2) = -ops[name].as<complex>();
        Hs_dn(s2, s1) = -conj(ops[name].as<complex>());
      }
    }
    arma::vec seigs_dn;
    arma::eig_sym(seigs_dn, Hs_dn);

    for (int nup = 0; nup <= nsites; ++nup)
      for (int ndn = 0; ndn <= nsites; ++ndn) {

        // Compute exact gs energy
        double e0_exact = 0;
        for (int i = 0; i < nup; ++i)
          e0_exact += seigs_up(i);
        for (int i = 0; i < ndn; ++i)
          e0_exact += seigs_dn(i);

        auto block3 = Electron(nsites, nup, ndn);
        auto H = matrixC(ops, block3, block3);
        REQUIRE(H.is_hermitian(1e-8));
        arma::vec evecs;
        arma::eig_sym(evecs, H);
        double e0 = evecs(0);
        // printf("nup: %d, ndn: %d, dim: %d, e0: %f, e0_exact: %f\n", nup, ndn,
        //        (int)evecs.size(), e0, e0_exact);
        REQUIRE(close(e0_exact, e0));
      }
  }

  // Test Heisenberg terms at half-filling
  for (int nsites = 2; nsites <= 6; ++nsites) {
    Log("electron_matrix: Heisenberg all-to-all comparison test, N={}",
        nsites);
    int nup = nsites / 2;
    int ndn = nsites - nup;
    auto block_spinhalf = Spinhalf(nsites, nup);
    auto block_electron = Electron(nsites, nup, ndn);

    auto ops = testcases::spinhalf::HB_alltoall(nsites);
    auto ops_U = ops;
    ops_U += "U" * Op("HubbardU");
    ops_U["U"] = 999999; // gap out doubly occupied sites
    auto H_spinhalf = matrix(ops, block_spinhalf, block_spinhalf);
    auto H_electron = matrix(ops_U, block_electron, block_electron);
    REQUIRE(H_spinhalf.is_hermitian(1e-8));
    REQUIRE(H_electron.is_hermitian(1e-8));

    arma::vec eigs_spinhalf;
    arma::eig_sym(eigs_spinhalf, H_spinhalf);

    arma::vec eigs_electron;
    arma::eig_sym(eigs_electron, H_electron);

    for (uint64_t idx = 0; idx < eigs_spinhalf.size(); ++idx)
      REQUIRE(std::abs(eigs_spinhalf(idx) - eigs_electron(idx)) < 1e-5);
  }

  /////////////////////////////////////////////////////
  // Henry's MATLAB code test (tests Heisenberg terms)
  {
    Log("electron_matrix: U-hopping-HB full spectrum of Henry's Matlab code");

    auto [ops, eigs_correct] = randomAlltoAll4NoU();

    // Compute full spectrum in xdiag
    nsites = 4;
    std::vector<double> all_eigs;
    for (int nup = 0; nup <= nsites; ++nup)
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block = Electron(nsites, nup, ndn);
        auto H = matrix(ops, block, block);
        REQUIRE(H.is_hermitian(1e-8));
        arma::vec eigs;
        arma::eig_sym(eigs, H);
        for (auto eig : eigs)
          all_eigs.push_back(eig);
      }
    std::sort(all_eigs.begin(), all_eigs.end());
    REQUIRE(close(arma::vec(all_eigs), eigs_correct));

    std::tie(ops, eigs_correct) = randomAlltoAll4();
    all_eigs.clear();
    for (int nup = 0; nup <= nsites; ++nup)
      for (int ndn = 0; ndn <= nsites; ++ndn) {

        auto block = Electron(nsites, nup, ndn);
        auto H = matrix(ops, block, block);
        REQUIRE(H.is_hermitian(1e-8));
        arma::vec eigs;
        arma::eig_sym(eigs, H);

        for (auto eig : eigs)
          all_eigs.push_back(eig);
      }
    std::sort(all_eigs.begin(), all_eigs.end());
    REQUIRE(close(arma::vec(all_eigs), eigs_correct));
  }

  for (int N = 3; N <= 6; ++N) {
    Log.out("electron_matrix: random all-to-all complex exchange test, N={}",
            N);

    auto ops = xdiag::testcases::tj::tj_alltoall_complex(N);

    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {
        auto block = Electron(N, nup, ndn);
        auto H = matrixC(ops, block, block);
        REQUIRE(H.is_hermitian(1e-8));
      }

    // Set hoppings to zero
    for (int s1 = 0; s1 < N; ++s1)
      for (int s2 = s1 + 1; s2 < N; ++s2) {
        std::stringstream ss;
        ss << "T" << s1 << "_" << s2;
        std::string name = ss.str();
        ops[name] = 0.;
      }
    auto ops_U = ops;
    ops_U += "U" * Op("HubbardU");
    ops_U["U"] = 1000;

    // Check whether eigenvalues agree with HB model
    for (int nup = 0; nup <= N; ++nup) {
      int ndn = N - nup;
      auto block1 = Electron(N, nup, ndn);
      auto block2 = Spinhalf(N, nup);
      auto H1 = matrixC(ops_U, block1, block1);
      auto H2 = matrixC(ops, block2, block2);
      arma::vec eigs1;
      arma::eig_sym(eigs1, H1);

      arma::vec eigs2;
      arma::eig_sym(eigs2, H2);

      // Log("eigs1(0): {}, eigs2(0): {}", eigs1(0), eigs2(0));
      // Log("eigs1.size(): {}, eigs2.size(): {}", eigs1.size(), eigs2.size());

      arma::vec eigs1_sub = eigs1.subvec(0, eigs2.size() - 1);
      REQUIRE(close(eigs1_sub, eigs2, 1e-6, 1e-6));
    }
  }

  for (int N = 3; N <= 4; ++N) {
    Log.out("electron_matrix: random all-to-all complex exchange test Np "
            "<-> NoNp, N={}",
            N);
    auto ops = xdiag::testcases::tj::tj_alltoall_complex(N);
    test_electron_np_no_np_matrix(N, ops);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
}
