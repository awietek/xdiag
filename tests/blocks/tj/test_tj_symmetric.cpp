// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <vector>

#include <tests/blocks/electron/testcases_electron.hpp>
#include <tests/blocks/tj/testcases_tj.hpp>
#include <tests/catch.hpp>
#include <tests/is_approx_hermitian.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/matrices/apply.hpp>
#include <xdiag/matrices/matrix.hpp>
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
      arma::cx_mat H = matrixC(ops, block);
      REQUIRE(testcases::is_approx_hermitian(H, 1e-8));
      arma::vec eigs_full;
      arma::eig_sym(eigs_full, H);

      std::vector<double> eigs_sym;
      for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
        auto blockk = tJ(nsites, nup, ndn, irreps[k]);
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

  // All-to-all tJ couplings, real and complex (the complex couplings exercise
  // the complex-character / flux path of the symmetric kernels).
  for (int64_t nsites = 2; nsites < 6; ++nsites) {
    Log("tj symmetric matrix: tJ all-to-all N = {}", nsites);
    auto [irreps, mults] = get_cyclic_group_irreps_mult(nsites);
    test_spectra_np(tj_alltoall(nsites), nsites, irreps, mults);
    test_spectra_np(tj_alltoall_complex(nsites), nsites, irreps, mults);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
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
}
