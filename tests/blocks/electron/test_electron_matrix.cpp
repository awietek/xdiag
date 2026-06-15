// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <tests/blocks/electron/testcases_electron.hpp>
#include <tests/blocks/random_opsum_matrix.hpp>
#include <tests/blocks/spinhalf/testcases_spinhalf.hpp>
#include <tests/blocks/tj/testcases_tj.hpp>
#include <tests/catch.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/matrices/matrix.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

// On-site identity m12{s,s} == m1{s} * m1{s} for the named diagonal operators.
static void test_onsite(std::string op1, std::string op12) {
  for (int nsites = 2; nsites < 5; ++nsites) {
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto b = Electron(nsites, nup, ndn);
        for (int s = 0; s < nsites; ++s) {
          arma::mat m1 = matrix(Op(op1, s), b);
          arma::mat m12 = matrix(Op(op12, {s, s}), b);
          REQUIRE(isapprox(m12, arma::mat(m1 * m1)));
        }
      }
    }
  }
}

// Full (number-non-conserving) spectrum must equal the union of all number
// sectors' spectra.
static void test_electron_np_no_np_matrix(int nsites, OpSum ops) try {
  auto block_full = Electron(nsites);
  auto H_full = matrixC(ops, block_full, block_full);
  REQUIRE(H_full.is_hermitian(1e-12));
  arma::Col<double> eigs_full;
  arma::eig_sym(eigs_full, H_full);

  std::vector<double> all_eigs;
  for (int nup = 0; nup <= nsites; ++nup) {
    for (int ndn = 0; ndn <= nsites; ++ndn) {
      auto block = Electron(nsites, nup, ndn);
      auto H = matrixC(ops, block, block);
      REQUIRE(H.is_hermitian(1e-12));
      arma::Col<double> eigs;
      arma::eig_sym(eigs, H);
      for (auto eig : eigs) {
        all_eigs.push_back(eig);
      }
    }
  }
  std::sort(all_eigs.begin(), all_eigs.end());
  REQUIRE(isapprox(arma::Col(all_eigs), eigs_full));
}
XDIAG_CATCH

TEST_CASE("electron_matrix", "[electron]") try {
  using namespace xdiag::testcases::electron;

  // Hubbard chain reference spectrum (Weisse & Fehske), N=4, nup=3, ndn=2.
  Log("electron_matrix: Hubbard Weisse & Fehske spectrum");
  {
    int nsites = 4;
    double t = 1.0;
    double U = 5.0;
    auto block = Electron(nsites, 3, 2);
    OpSum ops;
    for (int i = 0; i < nsites; ++i) {
      ops += "T" * Op("Hop", {i, (i + 1) % nsites});
    }
    ops += "U" * Op("HubbardU");
    ops["T"] = t;
    ops["U"] = U;
    auto H1 = matrix(ops, block, block);
    REQUIRE(H1.is_hermitian(1e-12));

    double tp = t, tm = -t, UU = U, U2 = 2 * U;
    arma::Mat<double> H1_correct = {
        {U2, tm, 0., 0., tp, 0., tm, 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0., tm, 0., 0., 0., 0., 0.},
        {tm, U2, tm, tm, 0., tp, 0., tm, 0., 0., 0., 0.,
         0., 0., 0., 0., 0., 0., 0., tm, 0., 0., 0., 0.},
        {0., tm, U2, 0., tm, 0., 0., 0., tm, 0., 0., 0.,
         0., 0., 0., 0., 0., 0., 0., 0., tm, 0., 0., 0.},
        {0., tm, 0., UU, tm, 0., 0., 0., 0., tm, 0., 0.,
         0., 0., 0., 0., 0., 0., 0., 0., 0., tm, 0., 0.},
        {tp, 0., tm, tm, UU, tm, 0., 0., 0., 0., tm, 0.,
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0., tm, 0.},
        {0., tp, 0., 0., tm, UU, 0., 0., 0., 0., 0., tm,
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., tm},
        {tm, 0., 0., 0., 0., 0., U2, tm, 0., 0., tp, 0.,
         tm, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., tm, 0., 0., 0., 0., tm, UU, tm, tm, 0., tp,
         0., tm, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., tm, 0., 0., 0., 0., tm, UU, 0., tm, 0.,
         0., 0., tm, 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., tm, 0., 0., 0., tm, 0., U2, tm, 0.,
         0., 0., 0., tm, 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., tm, 0., tp, 0., tm, tm, U2, tm,
         0., 0., 0., 0., tm, 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., tm, 0., tp, 0., 0., tm, UU,
         0., 0., 0., 0., 0., tm, 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., tm, 0., 0., 0., 0., 0.,
         UU, tm, 0., 0., tp, 0., tm, 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., tm, 0., 0., 0., 0.,
         tm, U2, tm, tm, 0., tp, 0., tm, 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., tm, 0., 0., 0.,
         0., tm, UU, 0., tm, 0., 0., 0., tm, 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., tm, 0., 0.,
         0., tm, 0., U2, tm, 0., 0., 0., 0., tm, 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., tm, 0.,
         tp, 0., tm, tm, UU, tm, 0., 0., 0., 0., tm, 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., tm,
         0., tp, 0., 0., tm, U2, 0., 0., 0., 0., 0., tm},
        {tm, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         tm, 0., 0., 0., 0., 0., UU, tm, 0., 0., tp, 0.},
        {0., tm, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0., tm, 0., 0., 0., 0., tm, UU, tm, tm, 0., tp},
        {0., 0., tm, 0., 0., 0., 0., 0., 0., 0., 0., 0.,
         0., 0., tm, 0., 0., 0., 0., tm, U2, 0., tm, 0.},
        {0., 0., 0., tm, 0., 0., 0., 0., 0., 0., 0., 0.,
         0., 0., 0., tm, 0., 0., 0., tm, 0., UU, tm, 0.},
        {0., 0., 0., 0., tm, 0., 0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., tm, 0., tp, 0., tm, tm, U2, tm},
        {0., 0., 0., 0., 0., tm, 0., 0., 0., 0., 0., 0.,
         0., 0., 0., 0., 0., tm, 0., tp, 0., 0., tm, U2}};
    REQUIRE(isapprox(H1, H1_correct));
  }

  // Two-site Hubbard exact ground state energy.
  {
    auto block2 = Electron(2, 1, 1);
    for (int i = 0; i < 20; ++i) {
      double U = 1.234 * i;
      OpSum ops;
      ops += "T" * Op("Hop", {0, 1});
      ops += "U" * Op("HubbardU");
      ops["T"] = 1.0;
      ops["U"] = U;
      double e0_exact = 0.5 * (U - std::sqrt(U * U + 16));
      auto H = matrix(ops, block2, block2);
      REQUIRE(H.is_hermitian(1e-8));
      arma::Col<double> eigs;
      arma::eig_sym(eigs, H);
      REQUIRE(isapprox(e0_exact, eigs(0)));
    }
  }

  // Free fermions, random all-to-all (real).
  for (int nsites = 3; nsites < 7; ++nsites) {
    Log("electron_matrix: free fermion random all-to-all (real), N={}", nsites);
    OpSum ops = freefermion_alltoall(nsites);

    arma::Mat<double> Hs(nsites, nsites, arma::fill::zeros);
    for (auto const &[cpl, mono] : ops.plain()) {
      REQUIRE(mono.size() == 1);
      int s1 = mono[0][0], s2 = mono[0][1];
      double c = cpl.scalar().as<double>();
      Hs(s1, s2) = -c;
      Hs(s2, s1) = -c;
    }
    arma::vec seigs;
    arma::eig_sym(seigs, Hs);

    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        double e0_exact = 0;
        for (int i = 0; i < nup; ++i) {
          e0_exact += seigs(i);
        }
        for (int i = 0; i < ndn; ++i) {
          e0_exact += seigs(i);
        }
        auto block = Electron(nsites, nup, ndn);
        auto Hr = matrix(ops, block, block);
        REQUIRE(Hr.is_hermitian(1e-8));
        arma::vec eigsr;
        arma::eig_sym(eigsr, Hr);
        auto Hc = matrixC(ops, block, block);
        REQUIRE(Hc.is_hermitian(1e-8));
        arma::vec eigsc;
        arma::eig_sym(eigsc, Hc);
        REQUIRE(isapprox(eigsr, eigsc));
        REQUIRE(isapprox(e0_exact, eigsr(0)));
      }
    }
  }

  // Free fermions, random all-to-all (complex, up/dn different).
  for (int nsites = 3; nsites < 7; ++nsites) {
    Log("electron_matrix: free fermion random all-to-all (cplx), N={}", nsites);
    OpSum ops = freefermion_alltoall_complex_updn(nsites);

    arma::cx_mat Hs_up(nsites, nsites, arma::fill::zeros);
    arma::cx_mat Hs_dn(nsites, nsites, arma::fill::zeros);
    for (auto const &[cpl, mono] : ops.plain()) {
      int s1 = mono[0][0], s2 = mono[0][1];
      complex c = cpl.scalar().as<complex>();
      std::string type = mono[0].type();
      // Hopup/Hopdn: symmetric (-c, -c); HopupAsym/HopdnAsym: antisym (-c, +c).
      if (type == "Hopup") {
        Hs_up(s1, s2) += -c;
        Hs_up(s2, s1) += -c;
      } else if (type == "HopupAsym") {
        Hs_up(s1, s2) += -c;
        Hs_up(s2, s1) += c;
      } else if (type == "Hopdn") {
        Hs_dn(s1, s2) += -c;
        Hs_dn(s2, s1) += -c;
      } else if (type == "HopdnAsym") {
        Hs_dn(s1, s2) += -c;
        Hs_dn(s2, s1) += c;
      }
    }
    arma::vec seigs_up, seigs_dn;
    arma::eig_sym(seigs_up, Hs_up);
    arma::eig_sym(seigs_dn, Hs_dn);

    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        double e0_exact = 0;
        for (int i = 0; i < nup; ++i) {
          e0_exact += seigs_up(i);
        }
        for (int i = 0; i < ndn; ++i) {
          e0_exact += seigs_dn(i);
        }
        auto block = Electron(nsites, nup, ndn);
        auto H = matrixC(ops, block, block);
        REQUIRE(H.is_hermitian(1e-8));
        arma::vec eigs;
        arma::eig_sym(eigs, H);
        REQUIRE(isapprox(e0_exact, eigs(0)));
      }
    }
  }

  // Heisenberg terms at half filling: electron (with large U) vs. Spinhalf.
  for (int nsites = 2; nsites <= 6; ++nsites) {
    Log("electron_matrix: Heisenberg all-to-all comparison, N={}", nsites);
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

    arma::vec eigs_spinhalf, eigs_electron;
    arma::eig_sym(eigs_spinhalf, H_spinhalf);
    arma::eig_sym(eigs_electron, H_electron);
    for (uint64_t idx = 0; idx < eigs_spinhalf.size(); ++idx) {
      REQUIRE(std::abs(eigs_spinhalf(idx) - eigs_electron(idx)) < 1e-5);
    }
  }

  // Full spectrum vs. reference (Henry's MATLAB), hopping + U + Heisenberg.
  {
    Log("electron_matrix: U-hopping-HB full spectrum (reference)");
    int nsites = 4;
    auto [ops, eigs_correct] = randomAlltoAll4NoU();
    std::vector<double> all_eigs;
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block = Electron(nsites, nup, ndn);
        auto H = matrix(ops, block, block);
        REQUIRE(H.is_hermitian(1e-8));
        arma::vec eigs;
        arma::eig_sym(eigs, H);
        for (auto eig : eigs) {
          all_eigs.push_back(eig);
        }
      }
    }
    std::sort(all_eigs.begin(), all_eigs.end());
    REQUIRE(isapprox(arma::vec(all_eigs), eigs_correct));

    auto [opsU, eigs_correctU] = randomAlltoAll4();
    all_eigs.clear();
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block = Electron(nsites, nup, ndn);
        auto H = matrix(opsU, block, block);
        REQUIRE(H.is_hermitian(1e-8));
        arma::vec eigs;
        arma::eig_sym(eigs, H);
        for (auto eig : eigs) {
          all_eigs.push_back(eig);
        }
      }
    }
    std::sort(all_eigs.begin(), all_eigs.end());
    REQUIRE(isapprox(arma::vec(all_eigs), eigs_correctU));
  }

  // Complex tJ-style all-to-all: hermiticity over all number sectors.
  for (int N = 3; N <= 6; ++N) {
    Log("electron_matrix: random all-to-all complex exchange, N={}", N);
    auto ops = xdiag::testcases::tj::tj_alltoall_complex(N);
    for (int nup = 0; nup <= N; ++nup) {
      for (int ndn = 0; ndn <= N - nup; ++ndn) {
        auto block = Electron(N, nup, ndn);
        auto H = matrixC(ops, block, block);
        REQUIRE(H.is_hermitian(1e-8));
      }
    }
  }

  for (int N = 3; N <= 4; ++N) {
    Log("electron_matrix: complex exchange Np <-> NoNp, N={}", N);
    auto ops = xdiag::testcases::tj::tj_alltoall_complex(N);
    test_electron_np_no_np_matrix(N, ops);
  }

  // On-site identity SzSz{s,s} == Sz{s}^2. (NtotNtot / NupdnNupdn forbid
  // overlapping sites in the new architecture, so they have no same-site form.)
  test_onsite("Sz", "SzSz");

  // Same-site composites that DO require distinct sites must be rejected.
  {
    auto b = Electron(2);
    REQUIRE_THROWS(matrix(Op("NtotNtot", {0, 0}), b));
    REQUIRE_THROWS(matrix(Op("NupNdn", {0, 0}), b));
  }

  for (int nsites = 2; nsites < 5; ++nsites) {
    auto b = Electron(nsites);
    for (int s = 0; s < nsites; ++s) {
      arma::mat cdagup = matrix(Op("Cdagup", s), b);
      arma::mat cup = matrix(Op("Cup", s), b);
      arma::mat cdagdn = matrix(Op("Cdagdn", s), b);
      arma::mat cdn = matrix(Op("Cdn", s), b);

      // Nupdn{s} = Nup{s} Ndn{s}
      REQUIRE(isapprox(arma::mat(matrix(Op("Nup", s), b) *
                                 matrix(Op("Ndn", s), b)),
                       matrix(Op("Nupdn", s), b)));

      // Exchange{s,s} = 1/2 (S+_s S-_s + S-_s S+_s)
      arma::mat exch = 0.5 * (cdagup * cdn * cdagdn * cup +
                              cdagdn * cup * cdagup * cdn);
      REQUIRE(isapprox(exch, matrix(Op("Exchange", {s, s}), b)));

      // SdotS{s,s} = Exchange{s,s} + Sz_s Sz_s
      arma::mat sz = matrix(Op("Sz", s), b);
      REQUIRE(isapprox(arma::mat(exch + sz * sz),
                       matrix(Op("SdotS", {s, s}), b)));

      for (int s2 = 0; s2 < nsites; ++s2) {
        if (s2 == s) {
          continue; // these composites require distinct sites
        }
        REQUIRE(isapprox(arma::mat(matrix(Op("Nup", s), b) *
                                   matrix(Op("Ndn", s2), b)),
                         matrix(Op("NupNdn", {s, s2}), b)));
        REQUIRE(isapprox(arma::mat(matrix(Op("Ndn", s), b) *
                                   matrix(Op("Nup", s2), b)),
                         matrix(Op("NdnNup", {s, s2}), b)));
        REQUIRE(isapprox(arma::mat(matrix(Op("Nup", s), b) *
                                   matrix(Op("Nup", s2), b)),
                         matrix(Op("NupNup", {s, s2}), b)));
        REQUIRE(isapprox(arma::mat(matrix(Op("Ndn", s), b) *
                                   matrix(Op("Ndn", s2), b)),
                         matrix(Op("NdnNdn", {s, s2}), b)));
      }
    }
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
}

// Randomized cross-check of the full operator pipeline against naive matrix
// products on the full Fock space (shared harness): expansion, protected-op
// expansion in products, normal ordering, fermi signs, and kernel dispatch must
// agree with multiplying/adding the individual single-Op matrices.
TEST_CASE("electronrandomopsum", "[electron]") try {
  for (int nsites = 2; nsites < 5; ++nsites) {
    Log("Electron random OpSum matrix test: N = {}", nsites);
    for (uint32_t seed = 0; seed < 5; ++seed) {
      testcases::test_random_opsum_matrix(Electron(nsites), seed);
    }
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
}
