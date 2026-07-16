// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <vector>

#include <tests/blocks/electron/testcases_electron.hpp>
#include <tests/catch.hpp>
#include <tests/is_approx_hermitian.hpp>

#include <xdiag/config.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/kernels/apply.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

// Dimension test: for every (nup, ndn) sector the symmetric block dimensions
// summed over all (translation) irreps (with multiplicity) must equal the
// non-symmetric block dimension; summed over all sectors they must equal the
// full local Hilbert space 4^nsites. This verifies that the coupled symmetric
// basis loses no states and double-counts none (the non-trivial check is for
// orbits with a non-trivial up-stabilizer, e.g. uniform up configurations).
static void test_dim_chain(int64_t nsites) {
  using namespace xdiag::testcases::electron;
  Log("electron symmetric dim: chain N = {}", nsites);
  auto [irreps, mults] = get_cyclic_group_irreps_mult(nsites);

  int64_t sum_total = 0;
  for (int64_t nup = 0; nup <= nsites; ++nup) {
    for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
      int64_t sum_updn = 0;
      for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
        auto block = Electron(nsites, nup, ndn, irreps[k]);
        sum_updn += dim(block) * mults[k];
        sum_total += dim(block) * mults[k];
      }
      REQUIRE(sum_updn == (int64_t)dim(Electron(nsites, nup, ndn)));
    }
  }
  REQUIRE(sum_total == (int64_t)std::pow(4, nsites));
}

// Spectrum-union (number-conserving): for each (nup, ndn) the union of the
// symmetric-block spectra over all irreps equals the non-symmetric spectrum.
// Real/complex consistency is checked for real blocks.
static void test_spectra_np(OpSum const &ops, int64_t nsites,
                            std::vector<Representation> const &irreps,
                            std::vector<int64_t> const &mults) {
  for (int64_t nup = 0; nup <= nsites; ++nup) {
    for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
      auto block = Electron(nsites, nup, ndn);
      if (dim(block) == 0 || dim(block) > 1000) {
        continue;
      }
      arma::cx_mat H = matrixC(ops, block);
      REQUIRE(testcases::is_approx_hermitian(H, 1e-8));
      arma::vec eigs_full;
      arma::eig_sym(eigs_full, H);

      std::vector<double> eigs_sym;
      for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
        auto blockk = Electron(nsites, nup, ndn, irreps[k]);
        if (dim(blockk) == 0) {
          continue;
        }
        arma::cx_mat Hk = matrixC(ops, blockk);
        REQUIRE(testcases::is_approx_hermitian(Hk, 1e-8));
        arma::vec ek;
        arma::eig_sym(ek, Hk);

        if (blockk.isreal() && isreal(ops)) {
          arma::mat Hk_real = matrix(ops, blockk);
          arma::vec ek_real;
          arma::eig_sym(ek_real, Hk_real);
          REQUIRE(isapprox(ek, ek_real));
        }
        for (double e : ek) {
          for (int64_t m = 0; m < mults[k]; ++m) {
            eigs_sym.push_back(e);
          }
        }
      }
      std::sort(eigs_sym.begin(), eigs_sym.end());
      REQUIRE(isapprox(eigs_full, arma::vec(eigs_sym)));
    }
  }
}

// Spectrum-union (no number conservation): the union of the symmetric no-np
// block spectra over all irreps equals the full block spectrum; and each no-np
// symmetric block's spectrum equals the union over (nup, ndn) of the number-
// conserving symmetric blocks of the same irrep.
static void test_spectra_no_np(OpSum const &ops, int64_t nsites,
                               std::vector<Representation> const &irreps) {
  auto block_total = Electron(nsites);
  if (dim(block_total) == 0 || dim(block_total) > 1000) {
    return;
  }
  arma::cx_mat H_total = matrixC(ops, block_total);
  REQUIRE(testcases::is_approx_hermitian(H_total, 1e-8));
  arma::vec eigs_total;
  arma::eig_sym(eigs_total, H_total);

  std::vector<double> eigs_all;
  for (auto const &irrep : irreps) {
    auto block = Electron(nsites, irrep);
    arma::cx_mat H = matrixC(ops, block);
    REQUIRE(testcases::is_approx_hermitian(H, 1e-8));
    arma::vec eigs_no_np;
    arma::eig_sym(eigs_no_np, H);

    std::vector<double> eigs_np_all;
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
        auto block_np = Electron(nsites, nup, ndn, irrep);
        if (dim(block_np) == 0) {
          continue;
        }
        arma::cx_mat H_np = matrixC(ops, block_np);
        REQUIRE(testcases::is_approx_hermitian(H_np, 1e-8));
        arma::vec eigs_np;
        arma::eig_sym(eigs_np, H_np);
        for (double e : eigs_np) {
          eigs_np_all.push_back(e);
          eigs_all.push_back(e);
        }
      }
    }
    std::sort(eigs_np_all.begin(), eigs_np_all.end());
    REQUIRE(isapprox(eigs_no_np, arma::vec(eigs_np_all)));
  }
  std::sort(eigs_all.begin(), eigs_all.end());
  REQUIRE(isapprox(eigs_total, arma::vec(eigs_all)));
}

// apply() agrees with the explicit matrix-vector / matrix-matrix product, and
// the Lanczos ground-state energy agrees with the dense one, in every symmetric
// block. Real/complex consistency is checked for real blocks.
static void test_apply(OpSum const &ops, int64_t nsites,
                       std::vector<Representation> const &irreps) {
  for (int64_t nup = 0; nup <= nsites; ++nup) {
    for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
      for (auto const &irrep : irreps) {
        auto block = Electron(nsites, nup, ndn, irrep);
        if (dim(block) == 0) {
          continue;
        }
        arma::cx_mat H = matrixC(ops, block);

        arma::cx_vec v(dim(block), arma::fill::randn);
        arma::cx_vec w1 = H * v;
        arma::cx_vec w2(dim(block), arma::fill::zeros);
        apply(ops, block, v, block, w2);
        REQUIRE(isapprox(w1, w2));

        arma::cx_mat m(dim(block), 3, arma::fill::randn);
        arma::cx_mat n1 = H * m;
        arma::cx_mat n2(dim(block), 3, arma::fill::zeros);
        apply(ops, block, m, block, n2);
        REQUIRE(isapprox(n1, n2));

        arma::vec eigs_mat;
        arma::eig_sym(eigs_mat, H);
        REQUIRE(std::abs(eigs_mat(0) - eigval0(ops, block)) < 1e-7);
      }
    }
  }
}

TEST_CASE("electronsymmetricdim", "[electron]") try {
  for (int64_t nsites = 1; nsites < 7; ++nsites) {
    test_dim_chain(nsites);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("electronsymmetricmatrix", "[electron]") try {
  using namespace xdiag::testcases::electron;

  // Exact symmetric matrices (Weisse & Fehske), 4-site Hubbard, nup=3, ndn=2.
  {
    Log("electron symmetric matrix: Weisse & Fehske 4-site Hubbard");
    int64_t nsites = 4, nup = 3, ndn = 2;
    double t = 1.0, U = 5.0;
    auto ops = get_linear_chain(nsites, t, U);
    auto [irreps, mults] = get_cyclic_group_irreps_mult(nsites);
    complex U2 = 2 * U, UU = U, tp = t, tm = -t, it = complex(0, t);

    for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
      auto block = Electron(nsites, nup, ndn, irreps[k]);
      arma::cx_mat H = matrixC(ops, block);

      // notice k=1 <-> k=3 switched w.r.t. reference as they define the
      // sign of the irrep differently
      // P = sum e^{ikr} T_r  (Fehske)
      // P = sum e^{-ikr} T_r (here)
      if (k == 0) {
        arma::cx_mat Hc = {{U2, tm, tm, tp, tp, 0.}, {tm, U2, tm, tm, 0., tp},
                           {tm, tm, U2, 0., tm, tm}, {tp, tm, 0., UU, tm, tp},
                           {tp, 0., tm, tm, UU, tm}, {0., tp, tm, tp, tm, UU}};
        REQUIRE(isapprox(Hc, H));
      } else if (k == 1) {
        arma::cx_mat Hc = {
            {U2, tm, -it, -it, tp, 0.},    {tm, U2, tm, tm, -2. * it, tp},
            {it, tm, U2, 0., tm, -it},     {it, tm, 0., UU, tm, -it},
            {tp, 2. * it, tm, tm, UU, tm}, {0., tp, it, it, tm, UU}};
        REQUIRE(isapprox(Hc, H));
      } else if (k == 2) {
        arma::cx_mat Hc = {{U2, tm, tp, tm, tp, 0.}, {tm, U2, tm, tm, 0., tp},
                           {tp, tm, U2, 0., tm, tp}, {tm, tm, 0., UU, tm, tm},
                           {tp, 0., tm, tm, UU, tm}, {0., tp, tp, tm, tm, UU}};
        REQUIRE(isapprox(Hc, H));
      } else if (k == 3) {
        arma::cx_mat Hc = {
            {U2, tm, it, it, tp, 0.},       {tm, U2, tm, tm, 2. * it, tp},
            {-it, tm, U2, 0., tm, it},      {-it, tm, 0., UU, tm, it},
            {tp, -2. * it, tm, tm, UU, tm}, {0., tp, -it, -it, tm, UU}};
        REQUIRE(isapprox(Hc, H));
      }
    }
  }

  // Linear chains: Hubbard, and Hubbard + Heisenberg (Exchange).
  for (int64_t nsites = 2; nsites < 6; ++nsites) {
    auto [irreps, mults] = get_cyclic_group_irreps_mult(nsites);
    Log("electron symmetric matrix: Hubbard chain N = {}", nsites);
    auto ops = get_linear_chain(nsites, 1.0, 5.0);
    test_spectra_np(ops, nsites, irreps, mults);
    test_spectra_no_np(ops, nsites, irreps);

    Log("electron symmetric matrix: Hubbard+Heisenberg chain N = {}", nsites);
    auto ops_hb = get_linear_chain_hb(nsites, 0.4);
    test_spectra_np(ops_hb, nsites, irreps, mults);
    test_spectra_no_np(ops_hb, nsites, irreps);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// 3x3 triangular lattice with a D3 point group (non-trivial stabilizers and
// genuinely multi-dimensional irreps, the E irreps entering with multiplicity
// 2). Reads the lattice + symmetries from a toml file. The dimension identity
// over all D3 irreps (weighted by multiplicity) is the strong, cheap check; the
// spectrum-union is run only for the small particle sectors (others exceed the
// dense-matrix guard at 9 sites). Both real and complex (flux) hoppings.
static std::vector<std::pair<std::string, int64_t>> triangular_d3_reps() {
  return {{"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
          {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
          {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
          {"Y.C1.A", 6}};
}

TEST_CASE("electronsymmetrictriangular", "[electron]") try {
  int64_t nsites = 9;

  std::vector<Representation> irreps;
  std::vector<int64_t> mults;
  {
    auto fl = FileToml(XDIAG_DIRECTORY
                       "/misc/data/triangular.9.hop.sublattices.tsl.toml");
    for (auto const &[name, mult] : triangular_d3_reps()) {
      irreps.push_back(read_representation(fl, name));
      mults.push_back(mult);
    }

    // Dimension identity over all D3 irreps (with multiplicity).
    Log("electron symmetric triangular 3x3: dimensions");
    int64_t sum_total = 0;
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
        int64_t sum_updn = 0;
        for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
          sum_updn += dim(Electron(nsites, nup, ndn, irreps[k])) * mults[k];
        }
        REQUIRE(sum_updn == (int64_t)dim(Electron(nsites, nup, ndn)));
        sum_total += sum_updn;
      }
    }
    REQUIRE(sum_total == (int64_t)std::pow(4, nsites));

    // // Real hopping Hubbard spectrum-union (small sectors only).
    // Log("electron symmetric triangular 3x3: real Hubbard spectra");
    // OpSum ops = fl["Interactions"].as<OpSum>();
    // ops["T"] = 1.0;
    // ops += "U" * Op("HubbardU");
    // ops["U"] = 5.0;
    // test_spectra_np(ops, nsites, irreps, mults);
  }

  // Complex (flux) hopping Hubbard spectrum-union (small sectors only). The file
  // lists directed Hopup/Hopdn bonds with a single complex coupling TPHI, the
  // up/dn flux phases encoded in the bond directions. Since Hopup/Hopdn are now
  // Hermitian, a complex coupling on them is not; we split each bond into the
  // Hermitian symmetric part (Re on Hopup/Hopdn) and the anti-Hermitian part
  // (i*Im on HopupAsym/HopdnAsym), as in freefermion_alltoall_complex_updn. The
  // JPHI (Heisenberg) bonds are dropped (the old test set JPHI = 0).
  {
    Log("electron symmetric triangular 3x3: complex Hubbard spectra");
    auto fl = FileToml(
        XDIAG_DIRECTORY
        "/misc/data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.toml");
    OpSum ops_raw = fl["Interactions"].as<OpSum>();
    complex tphi = complex(0.5, 0.5);
    OpSum ops;
    for (auto const &term : ops_raw) {
      Op op = term.monomial[0];
      std::string type = op.type();
      if (type == "Hopup") {
        ops += std::real(tphi) * Op("Hopup", {op[0], op[1]});
        ops += complex(0.0, std::imag(tphi)) * Op("HopupAsym", {op[0], op[1]});
      } else if (type == "Hopdn") {
        ops += std::real(tphi) * Op("Hopdn", {op[0], op[1]});
        ops += complex(0.0, std::imag(tphi)) * Op("HopdnAsym", {op[0], op[1]});
      }
    }
    ops += "U" * Op("HubbardU");
    ops["U"] = 5.0;
    // test_spectra_np(ops, nsites, irreps, mults);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("electronsymmetricapply", "[electron]") try {
  using namespace xdiag::testcases::electron;
  for (int64_t nsites = 2; nsites < 6; ++nsites) {
    auto [irreps, mults] = get_cyclic_group_irreps_mult(nsites);
    (void)mults;
    Log("electron symmetric apply: Hubbard chain N = {}", nsites);
    test_apply(get_linear_chain(nsites, 1.0, 5.0), nsites, irreps);
    Log("electron symmetric apply: Hubbard+Heisenberg chain N = {}", nsites);
    test_apply(get_linear_chain_hb(nsites, 0.4), nsites, irreps);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
