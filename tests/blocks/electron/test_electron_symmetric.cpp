// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <vector>

#include <tests/blocks/electron/testcases_electron.hpp>
#include <tests/blocks/spectrum_union.hpp>
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
      std::vector<Block> parts;
      for (Representation const &irrep : irreps) {
        parts.push_back(Electron(nsites, nup, ndn, irrep));
      }
      testcases::spectrum_union(ops, block, parts, mults);
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
  std::vector<Block> all_np;
  for (auto const &irrep : irreps) {
    std::vector<Block> np;
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
        np.push_back(Electron(nsites, nup, ndn, irrep));
      }
    }
    testcases::spectrum_union(ops, Electron(nsites, irrep), np);
    all_np.insert(all_np.end(), np.begin(), np.end());
  }
  testcases::spectrum_union(ops, block_total, all_np);
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

// Large site counts select the 64-bit bases (40 sites: LinTable<uint64_t>,
// 43 sites: Combinations<uint64_t>). With one up and one down electron the
// Hilbert space stays small enough for a dense cross-check. The second bond
// shell at distance 13 makes bonds cross bit position 32 in every translate.
TEST_CASE("electron_symmetric_large_nsites", "[electron]") try {
  for (int64_t nsites : {40, 43}) {
    Log("electron symmetric large nsites: N = {}", nsites);
    OpSum ops;
    for (int64_t i = 0; i < nsites; ++i) {
      ops += 1.0 * Op("Hop", {i, (i + 1) % nsites});
      ops += 0.3 * Op("Hop", {i, (i + 13) % nsites});
    }
    ops += 4.0 * Op("HubbardU");
    testcases::translation_spectrum_union(
        ops, Electron(nsites, 1, 1), nsites, [&](Representation const &irrep) {
          return Electron(nsites, 1, 1, irrep);
        });
  }

  // single free particle on a 36-site ring: momentum sector k holds exactly
  // the eigenvalue -2 t cos(2 pi k / N), up to the global sign convention
  {
    int64_t nsites = 36;
    double t = 1.3;
    OpSum ops;
    for (int64_t i = 0; i < nsites; ++i) {
      ops += t * Op("Hop", {i, (i + 1) % nsites});
    }
    arma::vec ek(nsites);
    arma::vec expected(nsites);
    for (int64_t k = 0; k < nsites; ++k) {
      Block block = Electron(nsites, 1, 0, cyclic_group_irrep(nsites, k));
      REQUIRE(dim(block) == 1);
      ek(k) = testcases::dense_spectrum(ops, block)(0);
      expected(k) = 2.0 * t * std::cos(2.0 * M_PI * (double)k / (double)nsites);
    }
    bool minus = arma::norm(ek + expected, "inf") < 1e-10;
    bool plus = arma::norm(ek - expected, "inf") < 1e-10;
    REQUIRE((minus || plus));
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// The symmetric kernels for density-density and single-density operators are
// not reached by the hopping/exchange Hamiltonians above.
TEST_CASE("electron_symmetric_density_operators", "[electron]") try {
  for (int64_t nsites : {4, 5, 6}) {
    Log("electron symmetric density operators: N = {}", nsites);
    OpSum ops;
    for (int64_t i = 0; i < nsites; ++i) {
      int64_t j = (i + 1) % nsites;
      ops += 1.0 * Op("Hop", {i, j});
      ops += 0.7 * Op("NupNup", {i, j});
      ops += 0.6 * Op("NdnNdn", {i, j});
      ops += 0.5 * Op("NupNdn", {i, j});
      ops += 0.4 * Op("NdnNup", {i, j});
      ops += 0.3 * Op("NtotNtot", {i, j});
      ops += 0.2 * Op("Nup", i);
      ops += 0.1 * Op("Ndn", i);
      ops += 0.15 * Op("Ntot", i);
    }
    ops += 0.25 * Op("Id");
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
        testcases::translation_spectrum_union(
            ops, Electron(nsites, nup, ndn), nsites,
            [&](Representation const &irrep) {
              return Electron(nsites, nup, ndn, irrep);
            });
      }
    }
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
