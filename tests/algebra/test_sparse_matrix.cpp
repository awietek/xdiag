// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algebra/sparse/apply.hpp>
#include <xdiag/algebra/sparse/coo_matrix.hpp>
#include <xdiag/algebra/sparse/csc_matrix.hpp>
#include <xdiag/algebra/sparse/csr_matrix.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/common.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/utils/xdiag_show.hpp>

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
    }
  }
}

TEST_CASE("sparse_matrix", "[algebra]") try {
  Log("Test sparse_matrix");

  for (int64_t nsites = 2; nsites < 5; ++nsites) {
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

  for (int N = 3; N <= 5; ++N) {
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

  for (int nsites = 2; nsites <= 5; ++nsites) {
    Log(" sparse_matrix (Electron): Heisenberg all-to-all comparison test, "
        "N={}",
        nsites);
    int nup = nsites / 2;
    int ndn = nsites - nup;
    auto block = Electron(nsites, nup, ndn);
    auto ops = testcases::spinhalf::HB_alltoall(nsites);
    test_sparse_matrix(ops, block);
  }

  for (int nsites = 3; nsites < 6; ++nsites) {
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

  for (int nsites = 3; nsites < 6; ++nsites) {
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
    for (int nsites = 3; nsites < 6; ++nsites) {
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

  for (int64_t nsites = 2; nsites < 6; ++nsites) {
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
    for (int nsites = 2; nsites <= 6; ++nsites) {
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
}
