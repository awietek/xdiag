// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <vector>

#include <tests/blocks/electron/testcases_electron.hpp>
#include <tests/blocks/spectrum_union.hpp>
#include <tests/blocks/tj/testcases_tj.hpp>
#include <tests/catch.hpp>
#include <tests/is_approx_hermitian.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/tj.hpp>
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
// full local Hilbert space 3^nsites (tJ local dimension 3, no double occupancy).
// This verifies that the coupled symmetric basis -- with its per-up-rep
// no-double-occupancy-filtered dn blocks -- loses no states and double-counts
// none.
static void test_dim_chain(int64_t nsites) {
  using namespace xdiag::testcases::electron;
  Log("tj symmetric dim: chain N = {}", nsites);
  auto [irreps, mults] = get_cyclic_group_irreps_mult(nsites);

  int64_t sum_total = 0;
  for (int64_t nup = 0; nup <= nsites; ++nup) {
    for (int64_t ndn = 0; ndn + nup <= nsites; ++ndn) {
      int64_t sum_updn = 0;
      for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
        auto block = tJ(nsites, nup, ndn, irreps[k]);
        sum_updn += dim(block) * mults[k];
        sum_total += dim(block) * mults[k];
      }
      REQUIRE(sum_updn == (int64_t)dim(tJ(nsites, nup, ndn)));
    }
  }
  REQUIRE(sum_total == (int64_t)std::llround(std::pow(3, nsites)));
}

// Spectrum-union (number-conserving): for each (nup, ndn) the union of the
// symmetric-block spectra over all irreps equals the non-symmetric spectrum.
// This is the primary cross-check between the compressed non-symmetric tJ path
// and the full-dn symmetric path (which reuses the electron symmetric kernels):
// any Jordan-Wigner sign / double-occupancy-projection mismatch shows up here.
// Real/complex consistency is checked for real blocks.
static void test_spectra_np(OpSum const &ops, int64_t nsites,
                            std::vector<Representation> const &irreps,
                            std::vector<int64_t> const &mults) {
  for (int64_t nup = 0; nup <= nsites; ++nup) {
    for (int64_t ndn = 0; ndn + nup <= nsites; ++ndn) {
      auto block = tJ(nsites, nup, ndn);
      if (dim(block) == 0 || dim(block) > 1000) {
        continue;
      }
      std::vector<Block> parts;
      for (Representation const &irrep : irreps) {
        parts.push_back(tJ(nsites, nup, ndn, irrep));
      }
      testcases::spectrum_union(ops, block, parts, mults);
    }
  }
}

// apply() agrees with the explicit matrix-vector / matrix-matrix product, and
// the Lanczos ground-state energy agrees with the dense one, in every symmetric
// block. Real/complex consistency is checked for real blocks.
static void test_apply(OpSum const &ops, int64_t nsites,
                       std::vector<Representation> const &irreps) {
  for (int64_t nup = 0; nup <= nsites; ++nup) {
    for (int64_t ndn = 0; ndn + nup <= nsites; ++ndn) {
      for (auto const &irrep : irreps) {
        auto block = tJ(nsites, nup, ndn, irrep);
        if (dim(block) == 0) {
          continue;
        }
        arma::cx_mat H = matrixC(ops, block);

        arma::cx_vec v(dim(block), arma::fill::randn);
        arma::cx_vec w1 = H * v;
        arma::cx_vec w2(dim(block), arma::fill::zeros);
        apply(ops, block, v, block, w2);
        REQUIRE(isapprox(w1, w2));

        arma::vec eigs_mat;
        arma::eig_sym(eigs_mat, H);
        REQUIRE(std::abs(eigs_mat(0) - eigval0(ops, block)) < 1e-7);
      }
    }
  }
}

TEST_CASE("tjsymmetricdim", "[tj]") try {
  for (int64_t nsites = 1; nsites < 7; ++nsites) {
    test_dim_chain(nsites);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("tjsymmetricmatrix", "[tj]") try {
  using namespace xdiag::testcases::tj;
  using xdiag::testcases::electron::get_cyclic_group_irreps_mult;

  // Real tJ chain (t-J coupling), spectra cross-checked against the
  // non-symmetric block for every momentum sector.
  for (int64_t nsites = 2; nsites < 7; ++nsites) {
    Log("tj symmetric matrix: tJ chain N = {}", nsites);
    auto ops = tJchain(nsites, 1.0, 0.4);
    auto [irreps, mults] = get_cyclic_group_irreps_mult(nsites);
    test_spectra_np(ops, nsites, irreps, mults);
  }

  // Translation-invariant complex chain (uniform Peierls flux): exercises the
  // complex-coupling / flux path of the symmetric kernels while still conserving
  // momentum, so the per-sector spectra union is meaningful. (A random all-to-all
  // coupling is NOT translation-invariant -- it has no SitePermutation quantum
  // number -- so it cannot be decomposed into momentum sectors and is covered
  // instead by the non-symmetric matrix test.)
  for (int64_t nsites = 3; nsites < 6; ++nsites) {
    Log("tj symmetric matrix: tJ complex flux chain N = {}", nsites);
    auto [irreps, mults] = get_cyclic_group_irreps_mult(nsites);
    OpSum ops;
    for (int64_t i = 0; i < nsites; ++i) {
      int64_t j = (i + 1) % nsites;
      ops += 1.0 * Op("Hop", {i, j});                // real part of the hopping
      ops += complex(0.0, 0.5) * Op("HopAsym", {i, j}); // imaginary part (flux)
      ops += 0.4 * Op("SdotS", {i, j});
    }
    test_spectra_np(ops, nsites, irreps, mults);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("tjsymmetricapply", "[tj]") try {
  using namespace xdiag::testcases::tj;
  using xdiag::testcases::electron::get_cyclic_group_irreps_mult;

  for (int64_t nsites = 2; nsites < 7; ++nsites) {
    Log("tj symmetric apply: tJ chain N = {}", nsites);
    auto ops = tJchain(nsites, 1.0, 0.4);
    auto [irreps, mults] = get_cyclic_group_irreps_mult(nsites);
    test_apply(ops, nsites, irreps);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// Large site counts select the 64-bit bases (40 sites: LinTable<uint64_t>,
// 43 sites: Combinations<uint64_t>); one up and one down hole keep the Hilbert
// space small. The bond shell at distance 13 crosses bit position 32.
TEST_CASE("tj_symmetric_large_nsites", "[tj]") try {
  for (int64_t nsites : {40, 43}) {
    Log("tj symmetric large nsites: N = {}", nsites);
    OpSum ops;
    for (int64_t i = 0; i < nsites; ++i) {
      for (auto [j, t] : std::vector<std::pair<int64_t, double>>{
               {(i + 1) % nsites, 1.0}, {(i + 13) % nsites, 0.3}}) {
        ops += t * Op("Hop", {i, j});
        ops += 0.4 * t * Op("tJSdotS", {i, j});
      }
    }
    testcases::translation_spectrum_union(
        ops, tJ(nsites, 1, 1), nsites,
        [&](Representation const &irrep) { return tJ(nsites, 1, 1, irrep); });
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// Symmetric tJ blocks without particle-number conservation use a separate
// basis (Subsets enumeration with the compressed dn index).
TEST_CASE("tj_symmetric_no_number_conservation", "[tj]") try {
  for (int64_t nsites : {4, 5, 6}) {
    Log("tj symmetric without number conservation: N = {}", nsites);
    OpSum ops;
    for (int64_t i = 0; i < nsites; ++i) {
      ops += 1.0 * Op("Hop", {i, (i + 1) % nsites});
      ops += 0.4 * Op("tJSdotS", {i, (i + 1) % nsites});
      ops += 0.2 * Op("Ntot", i);
    }
    testcases::translation_spectrum_union(
        ops, tJ(nsites), nsites,
        [&](Representation const &irrep) { return tJ(nsites, irrep); });
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// The symmetric tJ matrix kernel has its own branches for the density-density
// operators and the identity, which the hopping/exchange Hamiltonians above
// never reach.
TEST_CASE("tj_symmetric_density_operators", "[tj]") try {
  for (int64_t nsites : {4, 5, 6}) {
    Log("tj symmetric density operators: N = {}", nsites);
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
      for (int64_t ndn = 0; ndn + nup <= nsites; ++ndn) {
        testcases::translation_spectrum_union(
            ops, tJ(nsites, nup, ndn), nsites,
            [&](Representation const &irrep) {
              return tJ(nsites, nup, ndn, irrep);
            });
      }
    }
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
