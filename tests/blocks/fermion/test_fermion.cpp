// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <vector>

#include <tests/blocks/fermion/testcases_fermion.hpp>
#include <tests/blocks/random_opsum_matrix.hpp>
#include <tests/catch.hpp>
#include <tests/is_approx_hermitian.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/correlation_matrix.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

using namespace xdiag;

static OpSum spinless_fermi_chain(int64_t nsites, double t, double V, double mu,
                                  bool pbc = false) {

  auto ops = OpSum();
  int end = pbc ? nsites : nsites - 1;
  for (int i = 0; i < end; ++i) {
    ops += t * Op("Hop", {i, (i + 1) % nsites});
    ops += V * Op("NN", {i, (i + 1) % nsites});
  }
  ops += mu * Op("TotalN");
  return ops;
}

// Spinless chain with a complex (Peierls-flux) hopping, to exercise complex
// coefficients together with complex irrep characters under symmetry.
static OpSum spinless_fermi_chain_flux(int64_t nsites, double t, double tflux,
                                       double V, double mu) {
  auto ops = OpSum();
  for (int i = 0; i < nsites; ++i) {
    ops += t * Op("Hop", {i, (i + 1) % nsites});
    ops += std::complex(0.0, tflux) * Op("HopAsym", {i, (i + 1) % nsites});
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
      REQUIRE(testcases::is_approx_hermitian(Hr, 1e-8));
      arma::vec eigsr;
      arma::eig_sym(eigsr, Hr);
      auto Hc = matrixC(ops, block);
      REQUIRE(testcases::is_approx_hermitian(Hc, 1e-8));
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
      REQUIRE(testcases::is_approx_hermitian(Hc, 1e-8));
      arma::vec eigsc;
      arma::eig_sym(eigsc, Hc);
      double e0 = eigsc(0);
      // printf("number: %d,  e0: %f, e0_exact: %f\n", number, e0, e0_exact);
      REQUIRE(isapprox(e0_exact, e0));
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
}

TEST_CASE("fermionsymmetry", "[fermion]") try {

  // Test whether quantum numbers work
  for (int nsites = 2; nsites < 5; ++nsites) {
    Log("Spinless fermion chain symmetry test: N = {}", nsites);
    auto ops = spinless_fermi_chain(nsites, 1.0, 2.0, 0.5, true);
    auto block = Fermion(nsites);
    double e0 = eigval0(ops, block);
    // Log("e0: {:.6f}", e0);
    double e0nmin = 999.999999;

    for (int number = 0; number <= nsites; ++number) {
      auto blockn = Fermion(nsites, number);
      auto [e0n, psi0n] = eig0(ops, blockn);
      // Log("n: {}, e0n: {:.6f}", number, e0n);
      if (e0n < e0nmin) {
        e0nmin = e0n;
      }

      auto CdagCn = correlation_matrixC(psi0n, "Cdag", "C");

      // XDIAG_SHOW(CdagCn);

      double e0nkmin = 99999.999;
      int kmin = 0;
      int ndegen = 0; // number of k-sectors reaching the lowest energy
      for (int k = 0; k < nsites; ++k) {
        auto blocknk = Fermion(nsites, number, cyclic_group_irrep(nsites, k));
        // XDIAG_SHOW(blocknk);
        if (dim(blocknk) > 0) {
          auto [e0nk, psi0nk] = eig0(ops, blocknk);
          // Log("n: {}, k: {}, e0nk: {:.6f}", number, k, e0nk);
          if (e0nk < e0nkmin - 1e-8) {
            e0nkmin = e0nk;
            kmin = k;
            ndegen = 1;
          } else if (isapprox(e0nk, e0nkmin)) {
            ++ndegen;
          }
        }
      }
      auto blocknk = Fermion(nsites, number, cyclic_group_irrep(nsites, kmin));
      auto [e0nk, psi0nk] = eig0(ops, blocknk);
      auto CdagCnk = correlation_matrixC(psi0nk, "Cdag", "C");
      // XDIAG_SHOW(CdagCnk);

      REQUIRE(isapprox(e0n, e0nkmin));
      // <Cdag_i C_j> is only basis-independent for a non-degenerate ground
      // state. When the lowest energy is reached in several momentum sectors,
      // the non-symmetric block returns an arbitrary combination of degenerate
      // ground states, so the correlation matrices need not agree.
      if (ndegen == 1) {
        REQUIRE(isapprox(CdagCn, CdagCnk, 1e-6, 1e-6));
      }
    }
    REQUIRE(isapprox(e0, e0nmin));

    double e0kmin = 99999.999;
    for (int k = 0; k < nsites; ++k) {
      auto blockk = Fermion(nsites, cyclic_group_irrep(nsites, k));
      if (dim(blockk) > 0) {
        double e0k = eigval0(ops, blockk);
        // Log("k: {}, e0k: {:.6f}", k, e0k);
        if (e0k < e0kmin) {
          e0kmin = e0k;
        }
      }
    }
    REQUIRE(isapprox(e0, e0kmin));
    // Log("");
  }
} catch (xdiag::Error e) {
  error_trace(e);
}

// Strong check of the fermionic symmetry-adapted basis: for every particle
// sector the union of the symmetric-block spectra over all momenta must equal
// the full (non-symmetric) spectrum, and the symmetric block dimensions must
// sum to the full dimension. The dimension identity verifies that fermionic
// zero-norm states are excluded exactly (none lost, none double-counted), which
// is only non-trivial when orbits have a non-trivial stabilizer (e.g. N=4, 6).
// Both a real and a complex (flux) Hamiltonian are checked.
TEST_CASE("fermionsymmetryspectrum", "[fermion]") try {
  for (int nsites = 2; nsites < 7; ++nsites) {
    Log("Spinless fermion full-spectrum symmetry test: N = {}", nsites);
    std::vector<OpSum> hamiltonians = {
        spinless_fermi_chain(nsites, 1.0, 2.0, 0.5, true),
        spinless_fermi_chain_flux(nsites, 1.0, 0.7, 2.0, 0.5)};

    for (auto const &ops : hamiltonians) {
      for (int number = 0; number <= nsites; ++number) {
        auto block = Fermion(nsites, number);
        arma::cx_mat H = matrixC(ops, block);
        REQUIRE(testcases::is_approx_hermitian(H, 1e-8));
        arma::vec eigs_full;
        arma::eig_sym(eigs_full, H);

        std::vector<double> eigs_sym;
        int64_t dimsum = 0;
        for (int k = 0; k < nsites; ++k) {
          auto blockk = Fermion(nsites, number, cyclic_group_irrep(nsites, k));
          int64_t d = dim(blockk);
          dimsum += d;
          if (d > 0) {
            arma::cx_mat Hk = matrixC(ops, blockk);
            REQUIRE(testcases::is_approx_hermitian(Hk, 1e-8));
            arma::vec ek;
            arma::eig_sym(ek, Hk);
            for (double e : ek) {
              eigs_sym.push_back(e);
            }
          }
        }
        REQUIRE(dimsum == (int64_t)dim(block));
        std::sort(eigs_sym.begin(), eigs_sym.end());
        REQUIRE(isapprox(eigs_full, arma::vec(eigs_sym)));
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
}

// A global (site-free) operator inside a product must expand: TotalN -> sum_i
// N{i}, so that Op("TotalN") * Op("N", j) == sum_i Op("N", i) * Op("N", j).
TEST_CASE("fermiontotalnproduct", "[fermion]") try {
  for (int nsites = 2; nsites < 5; ++nsites) {
    for (int number = 0; number <= nsites; ++number) {
      auto block = Fermion(nsites, number);
      for (int j = 0; j < nsites; ++j) {
        OpSum lhs(Op("TotalN") * Op("N", j));
        OpSum rhs;
        for (int i = 0; i < nsites; ++i) {
          rhs += 1.0 * (Op("N", i) * Op("N", j));
        }
        REQUIRE(isapprox(matrix(lhs, block), matrix(rhs, block)));
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
}

// Randomized cross-check of the full operator pipeline against naive matrix
// products (shared harness). Uses the full Fock space so every elementary op
// is an endomorphism.
TEST_CASE("fermionrandomopsum", "[fermion]") try {
  for (int nsites = 2; nsites < 6; ++nsites) {
    Log("Fermion random OpSum matrix test: N = {}", nsites);
    for (uint32_t seed = 0; seed < 5; ++seed) {
      testcases::test_random_opsum_matrix(Fermion(nsites), seed);
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
}
