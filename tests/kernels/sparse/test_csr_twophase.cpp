// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <numeric>
#include <vector>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/kernels/sparse/coo_matrix.hpp>
#include <xdiag/kernels/sparse/csr_matrix.hpp>
#include <xdiag/kernels/sparse/sparse_matrix_types.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/utils/error.hpp>

#include "../../blocks/electron/testcases_electron.hpp"
#include "../../blocks/fermion/testcases_fermion.hpp"
#include "../../blocks/spinhalf/testcases_spinhalf.hpp"
#include "../../blocks/tj/testcases_tj.hpp"

using namespace xdiag;
using namespace arma;

// The two-phase CSR build (csr_matrix_nnz -> allocate -> csr_matrix_fill),
// exactly as the Julia wrapper drives it, must reproduce the one-shot
// csr_matrix (which is itself checked against the dense matrix elsewhere). We
// compare the dense forms and the nnz count.
template <typename idx_t, typename coeff_t>
static void check_twophase(OpSum const &ops, Block const &block, int i0) {
  std::vector<int64_t> counts = csr_matrix_nnz<coeff_t>(ops, block, block);
  int64_t nrows = size(block);
  REQUIRE((int64_t)counts.size() == nrows);
  int64_t nnz = std::accumulate(counts.begin(), counts.end(), (int64_t)0);

  arma::Col<idx_t> rowptr((arma::uword)(nrows + 1));
  arma::Col<idx_t> col((arma::uword)nnz);
  arma::Col<coeff_t> data((arma::uword)nnz);
  csr_matrix_fill<idx_t, coeff_t>(ops, block, block, counts, rowptr.memptr(),
                                  col.memptr(), data.memptr(), (idx_t)i0);

  CSRMatrix<idx_t, coeff_t> twophase{(idx_t)nrows, (idx_t)nrows, rowptr,
                                     col,          data,         (idx_t)i0,
                                     false};
  CSRMatrix<idx_t, coeff_t> concrete =
      csr_matrix<idx_t, coeff_t>(ops, block, block, (idx_t)i0);
  REQUIRE(nnz == (int64_t)concrete.data.n_elem);
  REQUIRE(norm(to_dense(twophase) - to_dense(concrete)) < 1e-12);
}

// COO two-phase (coo_matrix_nnz -> allocate -> coo_matrix_fill) must reproduce
// the one-shot coo_matrix.
template <typename idx_t, typename coeff_t>
static void check_twophase_coo(OpSum const &ops, Block const &block, int i0) {
  int64_t nnz = coo_matrix_nnz<coeff_t>(ops, block, block);
  int64_t nrows = size(block);

  arma::Col<idx_t> row((arma::uword)nnz), col((arma::uword)nnz);
  arma::Col<coeff_t> data((arma::uword)nnz);
  coo_matrix_fill<idx_t, coeff_t>(ops, block, block, nnz, row.memptr(),
                                  col.memptr(), data.memptr(), (idx_t)i0);

  COOMatrix<idx_t, coeff_t> twophase{(idx_t)nrows, (idx_t)nrows, row,
                                     col,          data,         (idx_t)i0,
                                     false};
  COOMatrix<idx_t, coeff_t> concrete =
      coo_matrix<idx_t, coeff_t>(ops, block, block, (idx_t)i0);
  REQUIRE(nnz == (int64_t)concrete.data.n_elem);
  REQUIRE(norm(to_dense(twophase) - to_dense(concrete)) < 1e-12);
}

// All idx/coeff combinations and both index bases; the real coeff path only
// where the sector is real (otherwise the double CSR would drop imaginary
// entries).
static void check_all(OpSum const &ops, Block const &block) {
  for (int i0 = 0; i0 < 2; ++i0) {
    check_twophase<int32_t, complex>(ops, block, i0);
    check_twophase<int64_t, complex>(ops, block, i0);
    check_twophase_coo<int32_t, complex>(ops, block, i0);
    check_twophase_coo<int64_t, complex>(ops, block, i0);
    if (isreal(ops) && isreal(block)) {
      check_twophase<int32_t, double>(ops, block, i0);
      check_twophase<int64_t, double>(ops, block, i0);
      check_twophase_coo<int32_t, double>(ops, block, i0);
      check_twophase_coo<int64_t, double>(ops, block, i0);
    }
  }
}

TEST_CASE("csr_twophase_spinhalf", "[kernels]") try {
  for (int64_t N = 2; N <= 5; ++N) {
    // all-to-all Heisenberg for the non-symmetric blocks
    OpSum ops = testcases::spinhalf::HB_alltoall(N);
    check_all(ops, Spinhalf(N));
    for (int64_t nup = 0; nup <= N; ++nup) {
      check_all(ops, Spinhalf(N, nup));
    }
    // translation-invariant NN Heisenberg for the momentum sectors
    OpSum ops_nn;
    for (int64_t i = 0; i < N; ++i) {
      ops_nn += Op("SdotS", {i, (i + 1) % N});
    }
    for (int64_t k = 0; k < N; ++k) {
      Representation irrep = cyclic_group_irrep(N, k);
      for (int64_t nup = 0; nup <= N; ++nup) {
        check_all(ops_nn, Spinhalf(N, nup, irrep));
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("csr_twophase_fermion", "[kernels]") try {
  for (int64_t N = 2; N <= 5; ++N) {
    OpSum ops = testcases::fermion::freefermion_alltoall(N);
    check_all(ops, Fermion(N));
    for (int64_t nf = 0; nf <= N; ++nf) {
      check_all(ops, Fermion(N, nf));
    }
    // translation-invariant NN hopping for the momentum sectors
    OpSum ops_nn;
    for (int64_t i = 0; i < N; ++i) {
      ops_nn += Op("Hop", {i, (i + 1) % N});
    }
    for (int64_t nf = 0; nf <= N; ++nf) {
      for (int64_t k = 0; k < N; ++k) {
        check_all(ops_nn, Fermion(N, nf, cyclic_group_irrep(N, k)));
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("csr_twophase_electron", "[kernels]") try {
  for (int64_t N = 2; N <= 4; ++N) {
    OpSum ops = testcases::electron::freefermion_alltoall(N);
    OpSum opsc = testcases::electron::freefermion_alltoall_complex_updn(N);
    for (int64_t nup = 0; nup <= N; ++nup) {
      for (int64_t ndn = 0; ndn <= N; ++ndn) {
        check_all(ops, Electron(N, nup, ndn));
        check_all(opsc, Electron(N, nup, ndn)); // complex ops
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("csr_twophase_tj", "[kernels]") try {
  for (int N = 2; N <= 4; ++N) {
    OpSum ops = testcases::tj::tj_alltoall_complex(N);
    for (int64_t nup = 0; nup <= N; ++nup) {
      for (int64_t ndn = 0; nup + ndn <= N; ++ndn) {
        check_all(ops, tJ(N, nup, ndn));
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("csr_twophase_boson", "[kernels]") try {
  for (int64_t N = 2; N <= 4; ++N) {
    OpSum ops;
    for (int64_t i = 0; i < N; ++i) {
      ops += Op("Hop", {i, (i + 1) % N});
    }
    ops += 2.0 * Op("HubbardU");
    ops += 0.5 * Op("TotalN");
    for (int64_t d = 2; d <= 3; ++d) {
      check_all(ops, Boson(N, d));
      for (int64_t nb = 0; nb <= (d - 1) * N; ++nb) {
        check_all(ops, Boson(N, d, nb));
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}
