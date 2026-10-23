// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/kernels/sparse/apply.hpp>
#include <xdiag/kernels/sparse/coo_matrix.hpp>
#include <xdiag/kernels/sparse/csc_matrix.hpp>
#include <xdiag/kernels/sparse/csr_matrix.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

#include "../blocks/electron/testcases_electron.hpp"
#include "../blocks/spinhalf/testcases_spinhalf.hpp"
#include "../blocks/tj/testcases_tj.hpp"

using namespace xdiag;
using namespace arma;

template <typename op_t>
void test_sparse_matrix(op_t const &ops, Block const &block) {
  for (int i0 = 0; i0 < 2; ++i0) {
    {
      auto m1 = matrix(ops, block);
      auto m2 = to_dense(coo_matrix(OpSum(ops), block, i0));
      auto m3 = to_dense(coo_matrix<int32_t, double>(OpSum(ops), block, i0));
      auto m4 = to_dense(coo_matrix<int64_t, double>(OpSum(ops), block, i0));
      REQUIRE(norm(m1 - m2) < 1e-12);
      REQUIRE(norm(m1 - m3) < 1e-12);
      REQUIRE(norm(m1 - m4) < 1e-12);
    }

    {
      auto m1 = matrix(ops, block);
      auto m2 = to_dense(csr_matrix(OpSum(ops), block, i0));
      auto m3 = to_dense(csr_matrix<int32_t, double>(OpSum(ops), block, i0));
      auto m4 = to_dense(csr_matrix<int64_t, double>(OpSum(ops), block, i0));
      REQUIRE(norm(m1 - m2) < 1e-12);
      REQUIRE(norm(m1 - m3) < 1e-12);
      REQUIRE(norm(m1 - m4) < 1e-12);

      // test apply
      auto v = arma::vec(size(block), arma::fill::randu);
      auto w1 = m1 * v;
      auto w2 = apply(csr_matrix(OpSum(ops), block, i0), v);
      REQUIRE(norm(w1 - w2) < 1e-12);

      auto V = arma::mat(size(block), 7, arma::fill::randu);
      auto W1 = m1 * V;
      auto W2 = apply(csr_matrix(OpSum(ops), block, i0), V);
      REQUIRE(norm(W1 - W2) < 1e-12);
    }

    {
      auto m1 = matrix(ops, block);
      auto m2 = to_dense(csc_matrix(OpSum(ops), block, i0));
      auto m3 = to_dense(csc_matrix<int32_t, double>(OpSum(ops), block, i0));
      auto m4 = to_dense(csc_matrix<int64_t, double>(OpSum(ops), block, i0));
      REQUIRE(norm(m1 - m2) < 1e-12);
      REQUIRE(norm(m1 - m3) < 1e-12);
      REQUIRE(norm(m1 - m4) < 1e-12);

      // test arma conversion
      auto X = csc_matrix(OpSum(ops), block, 0);
      auto row = arma::conv_to<arma::uvec>::from(X.row);
      auto colptr = arma::conv_to<arma::uvec>::from(X.colptr);
      auto A = arma::sp_mat(row, colptr, X.data, X.nrows, X.ncols);
      auto D = conv_to<arma::mat>::from(A);
      // XDIAG_SHOW(D);
      // XDIAG_SHOW(m1);
      // XDIAG_SHOW(row);
      // XDIAG_SHOW(colptr);
      // XDIAG_SHOW(X.data);
      // REQUIRE(norm(D - m1) < 1e-12);
    }
  }
}

template <typename op_t>
void test_sparse_matrixC(op_t const &ops, Block const &block) {
  for (int i0 = 0; i0 < 2; ++i0) {
    {
      auto m1 = matrixC(ops, block);
      auto m2 = to_dense(coo_matrixC(OpSum(ops), block, i0));
      auto m3 = to_dense(coo_matrix<int32_t, complex>(OpSum(ops), block, i0));
      auto m4 = to_dense(coo_matrix<int64_t, complex>(OpSum(ops), block, i0));
      REQUIRE(norm(m1 - m2) < 1e-12);
      REQUIRE(norm(m1 - m3) < 1e-12);
      REQUIRE(norm(m1 - m4) < 1e-12);
    }

    {
      auto m1 = matrixC(ops, block);
      auto m2 = to_dense(csr_matrixC(OpSum(ops), block, i0));
      auto m3 = to_dense(csr_matrix<int32_t, complex>(OpSum(ops), block, i0));
      auto m4 = to_dense(csr_matrix<int64_t, complex>(OpSum(ops), block, i0));
      REQUIRE(norm(m1 - m2) < 1e-12);
      REQUIRE(norm(m1 - m3) < 1e-12);
      REQUIRE(norm(m1 - m4) < 1e-12);

      // test apply
      arma::Col<complex> v = arma::cx_vec(size(block), arma::fill::randu);
      auto w1 = m1 * v;
      auto w2 = xdiag::apply(csr_matrixC(OpSum(ops), block, i0), v);
      REQUIRE(norm(w1 - w2) < 1e-12);

      auto V = arma::cx_mat(size(block), 7, arma::fill::randu);
      auto W1 = m1 * V;
      auto W2 = xdiag::apply(csr_matrixC(OpSum(ops), block, i0), V);
      REQUIRE(norm(W1 - W2) < 1e-12);
    }

    {
      auto m1 = matrixC(ops, block);
      auto m2 = to_dense(csc_matrixC(OpSum(ops), block, i0));
      auto m3 = to_dense(csc_matrix<int32_t, complex>(OpSum(ops), block, i0));
      auto m4 = to_dense(csc_matrix<int64_t, complex>(OpSum(ops), block, i0));
      REQUIRE(norm(m1 - m2) < 1e-12);
      REQUIRE(norm(m1 - m3) < 1e-12);
      REQUIRE(norm(m1 - m4) < 1e-12);

      // test arma conversion
      auto X = csc_matrixC(OpSum(ops), block, 0);
      auto row = arma::conv_to<arma::uvec>::from(X.row);
      auto colptr = arma::conv_to<arma::uvec>::from(X.colptr);
      auto A = arma::sp_cx_mat(row, colptr, X.data, X.nrows, X.ncols);
      auto D = conv_to<arma::cx_mat>::from(A);
      // REQUIRE(norm(D - m1) < 1e-12);
    }
  }
}

TEST_CASE("sparse_matrix", "[kernels]") try {
  Log("Test sparse_matrix");

  for (int64_t nsites = 2; nsites < 4; ++nsites) {
    Log.out(" sparse_matrix (Electron): single Ntot operator, nsites: {}",
            nsites);
    for (int i = 0; i < nsites; ++i) {
      auto op = Op("Ntot", i);
      auto block = Electron(nsites);
      test_sparse_matrix(op, block);
      test_sparse_matrixC(op, block);

      for (int64_t nup = 0; nup <= nsites; ++nup) {
        for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
          auto block = Electron(nsites, nup, ndn);
          test_sparse_matrix(op, block);
          test_sparse_matrixC(op, block);
        }
      }
    }
  }

  for (int N = 3; N <= 4; ++N) {
    Log.out(" sparse_matrix (Electron): random all-to-all complex exchange "
            "test, N={}",
            N);

    auto ops = xdiag::testcases::tj::tj_alltoall_complex(N);

    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {
        auto block = Electron(N, nup, ndn);
        test_sparse_matrixC(ops, block);
      }
  }

  for (int nsites = 2; nsites <= 4; ++nsites) {
    Log(" sparse_matrix (Electron): Heisenberg all-to-all comparison test, "
        "N={}",
        nsites);
    int nup = nsites / 2;
    int ndn = nsites - nup;
    auto block = Electron(nsites, nup, ndn);
    auto ops = testcases::spinhalf::HB_alltoall(nsites);
    test_sparse_matrix(ops, block);
  }

  for (int nsites = 3; nsites < 5; ++nsites) {
    Log(" sparse_matrix (Electron): free fermion random all-to-all test, "
        "(real), "
        "N={}",
        nsites);
    auto ops = testcases::electron::freefermion_alltoall(nsites);
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block = Electron(nsites, nup, ndn);
        test_sparse_matrix(ops, block);
      }
    }
  }

  for (int nsites = 3; nsites < 5; ++nsites) {
    Log(" sparse_matrix (Electron): free fermion random all-to-all test, "
        "(cplx), "
        "N={}",
        nsites);
    auto ops = testcases::electron::freefermion_alltoall_complex_updn(nsites);
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block = Electron(nsites, nup, ndn);
        test_sparse_matrixC(ops, block);
      }
    }
  }

  {
    Log(" sparse_matrix (Electron): Weisse & Fehske matrix");
    int64_t nsites = 4;
    int64_t nup = 3;
    int64_t ndn = 2;
    double t = 1.0;
    double U = 5.0;
    auto ops = testcases::electron::get_linear_chain(nsites, t, U);
    auto [irreps, multiplicities] =
        testcases::electron::get_cyclic_group_irreps_mult(nsites);
    for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
      if ((k == 0) || (k == 2)) {
        auto block = Electron(nsites, nup, ndn, irreps[k]);
        test_sparse_matrix(ops, block);
        test_sparse_matrixC(ops, block);
      } else {
        auto block = Electron(nsites, nup, ndn, irreps[k]);
        test_sparse_matrixC(ops, block);
      }
    }
  }

  {
    for (int nsites = 3; nsites < 5; ++nsites) {
      Log(" sparse_matrix (Electron): Hubbard chain, nsites (+Heisenberg "
          "terms): N={}",
          nsites);
      auto ops = testcases::electron::get_linear_chain_hb(nsites, 0.4);
      auto [irreps, multiplicities] =
          testcases::electron::get_cyclic_group_irreps_mult(nsites);
      for (int nup = 0; nup <= nsites; ++nup) {
        for (int ndn = 0; ndn <= nsites; ++ndn) {
          for (auto irrep : irreps) {
            auto block = Electron(nsites, nup, ndn, irrep);
            test_sparse_matrixC(ops, block);
          }
        }
      }
    }
  }

  {
    for (int nsites = 3; nsites < 5; ++nsites) {
      Log(" sparse_matrix (Electron): creation / annihilation: N={}", nsites);
      std::vector<std::string> op_strs = {"Cdagup", "Cdagdn", "Cup", "Cdn"};

      for (int nup = 1; nup < nsites; ++nup) {
        for (int ndn = 1; ndn < nsites; ++ndn) {
          auto block = Electron(nsites, nup, ndn);
          for (int i = 0; i < nsites; ++i) {
            for (auto op_i_str : op_strs) {
              auto op = Op(op_i_str, i);
              test_sparse_matrix(op, block);
              test_sparse_matrixC(op, block);
            }
          }
        }
      }
    }
  }

  for (int64_t nsites = 2; nsites < 5; ++nsites) {
    Log.out(" sparse_matrix (tJ): HB chain, symmetric spectra test, nsites: {}",
            nsites);
    OpSum ops;
    for (int64_t s = 0; s < nsites; ++s) {
      ops += Op("tJSdotS", {s, (s + 1) % nsites});
    }
    auto [irreps, multiplicities] =
        testcases::electron::get_cyclic_group_irreps_mult(nsites);

    for (int64_t nup = 1; nup <= nsites; ++nup) {
      for (int64_t ndn = 1; ndn <= nsites - nup; ++ndn) {

        auto block = tJ(nsites, nup, ndn);
        if (block.size() < 1000) {
          test_sparse_matrix(ops, block);
          test_sparse_matrixC(ops, block);

          for (auto irrep : irreps) {
            auto block = tJ(nsites, nup, ndn, irrep);
            if (isreal(block) && isreal(ops)) {
              test_sparse_matrix(ops, block);
            } else {
              test_sparse_matrixC(ops, block);
            }
          }
        }
      }
    }
  }

  // Log(" sparse_matrix (tJ) 3x3 triangular staggered flux, complex");
  // std::string lfile = XDIAG_DIRECTORY
  //     "/misc/data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.toml";

  // auto fl = FileToml(lfile);
  // auto ops = fl["Interactions"].as<OpSum>();
  // std::vector<double> etas{0.0, 0.1, 0.2, 0.3};
  // auto group = fl["Symmetries"].as<PermutationGroup>();

  // std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
  //     {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
  //     {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
  //     {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
  //     {"Y.C1.A", 6}};

  // std::vector<Representation> irreps;
  // for (auto [name, mult] : rep_name_mult) {
  //   irreps.push_back(read_representation(fl, name));
  // }

  // int nsites = 9;
  // for (auto eta : etas) {
  //   Log("  eta: {:.2f}", eta);
  //   ops["TPHI"] = complex(cos(eta * M_PI), sin(eta * M_PI));
  //   ops["JPHI"] = complex(cos(2 * eta * M_PI), sin(2 * eta * M_PI));
  //   for (auto irrep : irreps) {
  //     for (int64_t nup = 0; nup <= nsites; ++nup) {
  //       for (int64_t ndn = 0; ndn <= nsites - nup; ++ndn) {
  //         auto block = tJ(nsites, nup, ndn, irrep);
  //         test_sparse_matrixC(ops, block);
  //       }
  //     }
  //   }
  // }

  {
    for (int nsites = 2; nsites <= 5; ++nsites) {
      Log(" sparse_matrix (Spinhalf): Heisenberg all-to-all, N={}", nsites);
      for (int nup = 0; nup <= nsites; ++nup) {
        auto ops = testcases::spinhalf::HB_alltoall(nsites);
        auto block = Spinhalf(nsites, nup);
        test_sparse_matrix(ops, block);
        test_sparse_matrixC(ops, block);
      }
    }
  }

  // {
  //   Log(" sparse_matrix (Spinhalf): Triangular 3x3");
  //   std::string lfile = XDIAG_DIRECTORY
  //       "/misc/data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.toml";
  //   auto fl = FileToml(lfile);
  //   auto ops = fl["Interactions"].as<OpSum>();
  //   ops["Jz1"] = 1.00;
  //   ops["Jz2"] = 0.23;
  //   ops["Jx1"] = 0.76;
  //   ops["Jx2"] = 0.46;
  //   std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
  //       {"Gamma.D6.A1", 1}, {"Gamma.D6.A2", 1}, {"Gamma.D6.B1", 1},
  //       {"Gamma.D6.B2", 1}, {"Gamma.D6.E1", 2}, {"Gamma.D6.E2", 2},
  //       {"K.D3.A1", 2},     {"K.D3.A2", 2},     {"K.D3.E", 4},
  //       {"Y.D1.A", 6},      {"Y.D1.B", 6}};
  //   std::vector<Representation> irreps;
  //   std::vector<int64_t> multiplicities;
  //   for (auto [name, mult] : rep_name_mult) {
  //     irreps.push_back(read_representation(fl, name));
  //     multiplicities.push_back(mult);
  //   }
  //   int nsites = 9;
  //   for (int nup = 0; nup <= nsites; ++nup) {
  //     for (auto irrep : irreps) {
  //       auto block = Spinhalf(nsites, nup, irrep);
  //       test_sparse_matrixC(ops, block);
  //     }
  //   }
  // }
} catch (Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("sparse_matrix_32bit_wrappers", "[kernels]") try {
  Log("Test sparse_matrix 32-bit wrappers");

  auto block = Spinhalf(4, 2);
  auto ops = testcases::spinhalf::HB_alltoall(4);
  OpSum opsC = complex(0.0, 1.0) * Op("SdotS", {0, 1});
  opsC += Op("SdotS", {1, 2});

  arma::mat m = matrix(ops, block);
  arma::cx_mat mC = matrixC(opsC, block);

  for (int32_t i0 = 0; i0 < 2; ++i0) {
    REQUIRE(norm(m - to_dense(coo_matrix_32(ops, block, i0))) < 1e-12);
    REQUIRE(norm(m - to_dense(csr_matrix_32(ops, block, i0))) < 1e-12);
    REQUIRE(norm(m - to_dense(csc_matrix_32(ops, block, i0))) < 1e-12);

    REQUIRE(norm(mC - to_dense(coo_matrixC_32(opsC, block, i0))) < 1e-12);
    REQUIRE(norm(mC - to_dense(csr_matrixC_32(opsC, block, i0))) < 1e-12);
    REQUIRE(norm(mC - to_dense(csc_matrixC_32(opsC, block, i0))) < 1e-12);

    // the named wrappers carry the index base and hermiticity flag
    CSRMatrix<int32_t, double> A = csr_matrix_32(ops, block, i0);
    REQUIRE(A.i0 == i0);
    REQUIRE(A.nrows == (int32_t)block.size());
    REQUIRE(A.ncols == (int32_t)block.size());
    REQUIRE(ishermitian(A));
    REQUIRE(isreal(A));

    CSCMatrix<int32_t, complex> B = csc_matrixC_32(opsC, block, i0);
    REQUIRE(B.i0 == i0);
    REQUIRE(!isreal(B));
    COOMatrix<int32_t, complex> C = coo_matrixC_32(opsC, block, i0);
    REQUIRE(C.i0 == i0);
    REQUIRE(!isreal(C));
  }
} catch (Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("sparse_matrix_two_blocks", "[kernels]") try {
  Log("Test sparse_matrix rectangular (block_in != block_out)");

  // S+ maps the nup=1 sector into the nup=2 sector of a 4-site chain
  auto block_in = Spinhalf(4, 1);
  auto block_out = Spinhalf(4, 2);
  OpSum ops;
  ops += 1.0 * Op("S+", 0);
  ops += 0.5 * Op("S+", 2);
  OpSum opsC;
  opsC += complex(0.0, 1.0) * Op("S+", 0);
  opsC += 0.5 * Op("S+", 2);

  arma::mat m = matrix(ops, block_in, block_out);
  arma::cx_mat mC = matrixC(opsC, block_in, block_out);
  REQUIRE(m.n_rows == block_out.size());
  REQUIRE(m.n_cols == block_in.size());

  for (int64_t i0 = 0; i0 < 2; ++i0) {
    int32_t i032 = (int32_t)i0;
    // int64 named wrappers
    REQUIRE(norm(m - to_dense(coo_matrix(ops, block_in, block_out, i0))) <
            1e-12);
    REQUIRE(norm(m - to_dense(csr_matrix(ops, block_in, block_out, i0))) <
            1e-12);
    REQUIRE(norm(m - to_dense(csc_matrix(ops, block_in, block_out, i0))) <
            1e-12);
    REQUIRE(norm(mC - to_dense(coo_matrixC(opsC, block_in, block_out, i0))) <
            1e-12);
    REQUIRE(norm(mC - to_dense(csr_matrixC(opsC, block_in, block_out, i0))) <
            1e-12);
    REQUIRE(norm(mC - to_dense(csc_matrixC(opsC, block_in, block_out, i0))) <
            1e-12);

    // int32 named wrappers
    REQUIRE(norm(m - to_dense(coo_matrix_32(ops, block_in, block_out, i032))) <
            1e-12);
    REQUIRE(norm(m - to_dense(csr_matrix_32(ops, block_in, block_out, i032))) <
            1e-12);
    REQUIRE(norm(m - to_dense(csc_matrix_32(ops, block_in, block_out, i032))) <
            1e-12);
    REQUIRE(norm(mC - to_dense(coo_matrixC_32(opsC, block_in, block_out,
                                              i032))) < 1e-12);
    REQUIRE(norm(mC - to_dense(csr_matrixC_32(opsC, block_in, block_out,
                                              i032))) < 1e-12);
    REQUIRE(norm(mC - to_dense(csc_matrixC_32(opsC, block_in, block_out,
                                              i032))) < 1e-12);

    // templated versions
    REQUIRE(norm(m - to_dense(coo_matrix<int32_t, double>(ops, block_in,
                                                          block_out, i032))) <
            1e-12);
    REQUIRE(norm(m - to_dense(csr_matrix<int64_t, double>(ops, block_in,
                                                          block_out, i0))) <
            1e-12);
    REQUIRE(norm(mC - to_dense(csc_matrix<int32_t, complex>(
                          opsC, block_in, block_out, i032))) < 1e-12);

    // a rectangular matrix is never flagged hermitian
    CSRMatrix<int64_t, double> A = csr_matrix(ops, block_in, block_out, i0);
    REQUIRE(A.nrows == block_out.size());
    REQUIRE(A.ncols == block_in.size());
    REQUIRE(!ishermitian(A));

    // apply maps vectors of size(block_in) to vectors of size(block_out)
    arma::vec v(block_in.size(), arma::fill::randu);
    REQUIRE(norm(m * v - xdiag::apply(A, v)) < 1e-12);
    arma::cx_vec vc(block_in.size(), arma::fill::randu);
    REQUIRE(norm(m * vc - xdiag::apply(A, vc)) < 1e-12);
    arma::mat V(block_in.size(), 3, arma::fill::randu);
    REQUIRE(norm(m * V - xdiag::apply(A, V)) < 1e-12);
    arma::cx_mat VC(block_in.size(), 3, arma::fill::randu);
    REQUIRE(norm(m * VC - xdiag::apply(A, VC)) < 1e-12);
  }
} catch (Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("sparse_matrix_apply_preallocated", "[kernels]") try {
  Log("Test sparse_matrix apply into preallocated output");

  auto block = Spinhalf(5, 2);
  int64_t D = block.size();
  auto ops = testcases::spinhalf::HB_alltoall(5);
  OpSum opsC = complex(0.3, 0.7) * Op("SdotS", {0, 1});
  opsC += Op("SdotS", {1, 3});
  opsC += Op("Sz", 2);

  arma::mat m = matrix(ops, block);
  arma::cx_mat mC = matrixC(opsC, block);

  arma::vec v(D, arma::fill::randu);
  arma::cx_vec vc(D, arma::fill::randu);
  arma::mat V(D, 4, arma::fill::randu);
  arma::cx_mat VC(D, 4, arma::fill::randu);

  for (int64_t i0 = 0; i0 < 2; ++i0) {
    int32_t i032 = (int32_t)i0;
    auto A64 = csr_matrix(ops, block, i0);
    auto A32 = csr_matrix_32(ops, block, i032);
    auto C64 = csr_matrixC(opsC, block, i0);
    auto C32 = csr_matrixC_32(opsC, block, i032);

    // real matrix, real input
    {
      arma::vec w(D);
      xdiag::apply(A64, v, w);
      REQUIRE(norm(m * v - w) < 1e-12);
      xdiag::apply(A32, v, w);
      REQUIRE(norm(m * v - w) < 1e-12);
      arma::mat W(D, 4);
      xdiag::apply(A64, V, W);
      REQUIRE(norm(m * V - W) < 1e-12);
      xdiag::apply(A32, V, W);
      REQUIRE(norm(m * V - W) < 1e-12);
    }

    // real matrix, complex input (mixed precision)
    {
      arma::cx_vec w(D);
      xdiag::apply(A64, vc, w);
      REQUIRE(norm(m * vc - w) < 1e-12);
      xdiag::apply(A32, vc, w);
      REQUIRE(norm(m * vc - w) < 1e-12);
      arma::cx_mat W(D, 4);
      xdiag::apply(A64, VC, W);
      REQUIRE(norm(m * VC - W) < 1e-12);
      xdiag::apply(A32, VC, W);
      REQUIRE(norm(m * VC - W) < 1e-12);

      // allocating overloads
      REQUIRE(norm(m * vc - xdiag::apply(A32, vc)) < 1e-12);
      REQUIRE(norm(m * VC - xdiag::apply(A32, VC)) < 1e-12);
    }

    // complex matrix, complex input
    {
      arma::cx_vec w(D);
      xdiag::apply(C64, vc, w);
      REQUIRE(norm(mC * vc - w) < 1e-12);
      xdiag::apply(C32, vc, w);
      REQUIRE(norm(mC * vc - w) < 1e-12);
      arma::cx_mat W(D, 4);
      xdiag::apply(C64, VC, W);
      REQUIRE(norm(mC * VC - W) < 1e-12);
      xdiag::apply(C32, VC, W);
      REQUIRE(norm(mC * VC - W) < 1e-12);

      REQUIRE(norm(mC * vc - xdiag::apply(C32, vc)) < 1e-12);
      REQUIRE(norm(mC * VC - xdiag::apply(C32, VC)) < 1e-12);
    }

    // the kernels:: dispatch helper refuses real vectors for complex matrices
    {
      arma::vec w(D);
      REQUIRE_THROWS_AS(kernels::apply(C64, v, w), Error);
      arma::cx_vec wc(D);
      kernels::apply(C64, vc, wc);
      REQUIRE(norm(mC * vc - wc) < 1e-12);
      kernels::apply(A64, v, w);
      REQUIRE(norm(m * v - w) < 1e-12);
    }
  }
} catch (Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("sparse_matrix_apply_errors", "[kernels]") try {
  Log("Test sparse_matrix apply error handling");

  auto block = Spinhalf(4, 2);
  int64_t D = block.size();
  auto ops = testcases::spinhalf::HB_alltoall(4);
  auto A = csr_matrix(ops, block);

  // mismatched vector / matrix dimensions
  {
    arma::vec v(D + 1, arma::fill::randu);
    REQUIRE_THROWS_AS(xdiag::apply(A, v), Error);
    arma::vec vin(D, arma::fill::randu);
    arma::vec wout(D - 1);
    REQUIRE_THROWS_AS(xdiag::apply(A, vin, wout), Error);

    arma::mat V(D - 1, 2, arma::fill::randu);
    REQUIRE_THROWS_AS(xdiag::apply(A, V), Error);
    arma::mat Vin(D, 2, arma::fill::randu);
    arma::mat Wout(D, 3);
    REQUIRE_THROWS_AS(xdiag::apply(A, Vin, Wout), Error);

    arma::cx_vec vc(D + 1, arma::fill::randu);
    REQUIRE_THROWS_AS(xdiag::apply(A, vc), Error);
    arma::cx_mat VC(D + 1, 2, arma::fill::randu);
    REQUIRE_THROWS_AS(xdiag::apply(A, VC), Error);
  }

  // an invalid index base is rejected
  {
    auto B = A;
    B.i0 = 2;
    arma::vec v(D, arma::fill::randu);
    REQUIRE_THROWS_AS(xdiag::apply(B, v), Error);
    arma::mat V(D, 2, arma::fill::randu);
    REQUIRE_THROWS_AS(xdiag::apply(B, V), Error);
  }

  // an invalid index base when building is rejected
  {
    REQUIRE_THROWS_AS(csr_matrix(ops, block, 2), Error);
    REQUIRE_THROWS_AS(csc_matrix(ops, block, -1), Error);
    REQUIRE_THROWS_AS(coo_matrix(ops, block, 2), Error);
  }
} catch (Error const &e) {
  error_trace(e);
  throw;
}
