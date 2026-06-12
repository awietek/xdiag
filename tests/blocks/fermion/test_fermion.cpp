// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/blocks/fermion/testcases_fermion.hpp>
#include <tests/catch.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/matrices/matrix.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/correlation_matrix.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

using namespace xdiag;

static OpSum spinless_fermi_chain(int64_t nsites, double t, double V,
                                  double mu) {

  auto ops = OpSum();
  for (int i = 0; i < nsites - 1; ++i) {
    ops += t * Op("Hop", {i, (i + 1) % nsites});
    ops += V * Op("NN", {i, (i + 1) % nsites});
  }
  ops += mu * Op("TotalN");
  return ops;
}

TEST_CASE("fermion", "[fermion]") try {
  Log("Fermion DMRG comparison");
  {
    Fermion block;
    REQUIRE_THROWS(block = Fermion(-1));
    REQUIRE_THROWS(block = Fermion(2, 3));
    REQUIRE_THROWS(block = Fermion(2, -1));
  }

  int64_t nsites = 6;
  double t = 1.0;
  auto block = Fermion(nsites);

  Log("Spinless Fermion chain");
  {
    double mu = 0.5;
    double V = 2.0;
    auto ops = spinless_fermi_chain(nsites, t, V, mu);
    auto [e0, psi0] = eig0(ops, block);
    double e0_dmrg = -1.9028099909816971;
    Log("e0: {:.16f}, e0_dmrg: {:.16f}", e0, e0_dmrg);
    REQUIRE(isapprox(e0, e0_dmrg, 1e-8, 1e-8));

    arma::mat CdagC_dmrg = {
        {0.277133, 0.356898, 0.223408, 0.0378285, -0.0846251, -0.101831},
        {0.356898, 0.472122, 0.314252, 0.0905766, -0.0485061, -0.0846251},
        {0.223408, 0.314252, 0.250745, 0.150337, 0.0905766, 0.0378285},
        {0.0378285, 0.0905766, 0.150337, 0.250745, 0.314252, 0.223408},
        {-0.0846251, -0.0485061, 0.0905766, 0.314252, 0.472122, 0.356898},
        {-0.101831, -0.0846251, 0.0378285, 0.223408, 0.356898, 0.277133}};
    auto CdagC = correlation_matrix(psi0, "Cdag", "C");
    REQUIRE(isapprox(CdagC_dmrg, CdagC, 1e-5, 1e-5));
  }

  {
    double mu = -0.2;
    double V = 3.7;
    auto ops = spinless_fermi_chain(nsites, t, V, mu);
    auto [e0, psi0] = eig0(ops, block);
    double e0_dmrg = -3.252522016281353;
    Log("e0: {:.16f}, e0_dmrg: {:.16f}", e0, e0_dmrg);
    REQUIRE(isapprox(e0, e0_dmrg, 1e-8, 1e-8));

    arma::mat CdagC_dmrg = {
        {0.294418, 0.369971, 0.219764, 0.0388638, -0.0693993, -0.0916508},
        {0.369971, 0.479203, 0.300678, 0.0809754, -0.0323916, -0.0693993},
        {0.219764, 0.300678, 0.226379, 0.124042, 0.0809754, 0.0388638},
        {0.0388638, 0.0809754, 0.124042, 0.226379, 0.300678, 0.219764},
        {-0.0693993, -0.0323916, 0.0809754, 0.300678, 0.479203, 0.369971},
        {-0.0916508, -0.0693993, 0.0388638, 0.219764, 0.369971, 0.294418}};
    auto CdagC = correlation_matrix(psi0, "Cdag", "C");
    REQUIRE(isapprox(CdagC_dmrg, CdagC, 1e-5, 1e-5));
  }
} catch (xdiag::Error e) {
  error_trace(e);
}

TEST_CASE("fermionanticommutation", "[fermion]") try {
  Log("Fermion anti-commutation relations");
  int nsites = 4;
  auto b = Fermion(nsites);
  int dim = b.dim();

  for (int i = 0; i < nsites; ++i) {
    for (int j = 0; j < nsites; ++j) {
      arma::mat cdagi = matrix(Op("Cdag", i), b);
      arma::mat cdagj = matrix(Op("Cdag", j), b);
      arma::mat ci = matrix(Op("C", i), b);
      arma::mat cj = matrix(Op("C", j), b);
      arma::mat id = arma::eye(dim, dim);
      arma::mat zeros = arma::mat(dim, dim, arma::fill::zeros);
      arma::mat cdagcdag = cdagi * cdagj + cdagj * cdagi;
      arma::mat cdagc = cdagi * cj + cj * cdagi;
      arma::mat ccdag = ci * cdagj + cdagj * ci;
      arma::mat cc = ci * cj + cj * ci;

      if (i == j) {
        REQUIRE(isapprox(cdagc, id));
        REQUIRE(isapprox(ccdag, id));
        REQUIRE(isapprox(cdagcdag, zeros));
        REQUIRE(isapprox(cc, zeros));
      } else {
        REQUIRE(isapprox(cdagc, zeros));
        REQUIRE(isapprox(ccdag, zeros));
        REQUIRE(isapprox(cdagcdag, zeros));
        REQUIRE(isapprox(cc, zeros));
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
}

TEST_CASE("randomfreefermionsreal", "[fermion]") try {
  for (int nsites = 3; nsites < 7; ++nsites) {
    Log("electron_matrix: free fermion random all-to-all test, (real), N={}",
        nsites);

    OpSum ops = testcases::fermion::freefermion_alltoall(nsites);

    // Create single particle matrix
    arma::Mat<double> Hs(nsites, nsites, arma::fill::zeros);
    for (auto [cpl, mono] : ops) {
      REQUIRE(mono.size() == 1);
      int s1 = mono[0][0];
      int s2 = mono[0][1];
      double c = cpl.scalar().as<double>();
      Hs(s1, s2) = -c;
      Hs(s2, s1) = -c;
    }

    arma::vec seigs;
    arma::eig_sym(seigs, Hs);
    for (int number = 0; number <= nsites; ++number) {

      // Compute exact gs energy
      double e0_exact = 0;
      for (int i = 0; i < number; ++i) {
        e0_exact += seigs(i);
      }

      auto block = Fermion(nsites, number);
      auto Hr = matrix(ops, block);
      REQUIRE(Hr.is_hermitian(1e-8));
      arma::vec eigsr;
      arma::eig_sym(eigsr, Hr);
      auto Hc = matrixC(ops, block);
      REQUIRE(Hc.is_hermitian(1e-8));
      arma::vec eigsc;
      arma::eig_sym(eigsc, Hc);
      REQUIRE(isapprox(eigsr, eigsc));
      double e0 = eigsr(0);
      REQUIRE(isapprox(e0_exact, e0));
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
}

TEST_CASE("randomfreefermionscomplex", "[fermion]") try {
  for (int nsites = 3; nsites < 7; ++nsites) {
    Log("electron_matrix: free fermion random all-to-all test, (cplx), N={}",
        nsites);

    OpSum ops = testcases::fermion::freefermion_alltoall_complex(nsites);

    // Create single particle matrix
    arma::cx_mat Hs(nsites, nsites, arma::fill::zeros);
    for (auto [cpl, mono] : ops) {
      REQUIRE(mono.size() == 1);
      int s1 = mono[0][0];
      int s2 = mono[0][1];
      complex c = cpl.scalar().as<complex>();
      Hs(s1, s2) += -c;
      Hs(s2, s1) += -std::conj(c);
    }

    arma::vec seigs;
    arma::eig_sym(seigs, Hs);
    for (int number = 0; number <= nsites; ++number) {

      // Compute exact gs energy
      double e0_exact = 0;
      for (int i = 0; i < number; ++i) {
        e0_exact += seigs(i);
      }

      auto block = Fermion(nsites, number);
      auto Hc = matrixC(ops, block);
      REQUIRE(Hc.is_hermitian(1e-8));
      arma::vec eigsc;
      arma::eig_sym(eigsc, Hc);
      double e0 = eigsc(0);
      printf("number: %d,  e0: %f, e0_exact: %f\n", number, e0, e0_exact);
      REQUIRE(isapprox(e0_exact, e0));
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
}
