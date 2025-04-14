// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include "testcases_spinhalf.hpp"
#include <iostream>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/utils/xdiag_show.hpp>

using namespace xdiag;

static void test_onsite(std::string op1, std::string op12) {
  for (int nsites = 2; nsites < 5; ++nsites) {
    for (int nup = 0; nup <= nsites; ++nup) {
      auto b = Spinhalf(nsites, nup);
      for (int s = 0; s < nsites; ++s) {
        arma::mat m1 = matrix(Op(op1, s), b);
        arma::mat m12 = matrix(Op(op12, {s, s}), b);
        REQUIRE(isapprox(m12, m1 * m1));
      }
    }
  }
}

TEST_CASE("spinhalf_matrix", "[spinhalf]") try {
  using namespace xdiag::testcases::spinhalf;

  {
    Log("spinhalf_matrix: Heisenberg chain test, J=1.0, N=2,..,6");
    for (int nsites = 2; nsites <= 6; ++nsites) {
      for (int nup = 0; nup <= nsites; ++nup) {
        auto [ops, exact_eigs] = HBchain_fullspectrum_nup(nsites, nup);
        auto block = Spinhalf(nsites, nup);
        auto H = matrix(ops, block, block);
        REQUIRE(H.is_hermitian(1e-7));
        arma::vec eigs;
        arma::eig_sym(eigs, H);
        // XDIAG_SHOW(ops);
        // H.print();
        // exact_eigs.print();
        REQUIRE(isapprox(eigs, exact_eigs));
      }
    }
  }

  {
    Log("spinhalf_matrix: Heisenberg all-to-all tJ comparison");
    for (int nsites = 2; nsites <= 6; ++nsites)
      for (int nup = 0; nup <= nsites; ++nup) {
        auto ops = HB_alltoall(nsites);
        auto block = Spinhalf(nsites, nup);
        auto block_tJ = tJ(nsites, nup, nsites - nup);
        auto H = matrix(ops, block, block);
        auto H_tJ = matrix(ops, block_tJ, block_tJ);
        REQUIRE(H.is_hermitian());
        REQUIRE(H_tJ.is_hermitian());

        arma::vec eigs;
        arma::eig_sym(eigs, H);

        arma::vec eigs_tJ;
        arma::eig_sym(eigs_tJ, H_tJ);
        REQUIRE(isapprox(eigs, eigs_tJ));
      }
  }

  {
    Log("spinhalf_matrix: Heisenberg all-to-all Sz <-> NoSz comparison");
    for (int nsites = 2; nsites <= 6; ++nsites) {
      auto ops = HB_alltoall(nsites);
      auto block_no_sz = Spinhalf(nsites);
      auto H_no_sz = matrix(ops, block_no_sz, block_no_sz);
      REQUIRE(H_no_sz.is_hermitian(1e-8));
      arma::vec eigs_no_sz;
      arma::eig_sym(eigs_no_sz, H_no_sz);

      std::vector<double> eigs_sz_all;

      for (int nup = 0; nup <= nsites; ++nup) {
        auto block_sz = Spinhalf(nsites, nup);
        auto H_sz = matrix(ops, block_sz, block_sz);
        REQUIRE(H_sz.is_hermitian(1e-7));
        arma::vec eigs_sz;
        arma::eig_sym(eigs_sz, H_sz);

        for (auto eig : eigs_sz)
          eigs_sz_all.push_back(eig);
      }
      std::sort(eigs_sz_all.begin(), eigs_sz_all.end());

      REQUIRE(isapprox(eigs_no_sz, arma::vec(eigs_sz_all)));
    }
  }

  {
    Log("spinhalf_matrix: triangular N=12 complex exchange");
    int nsites = 12;
    int nup = 6;
    std::vector<double> etas = {0.00, 0.01, 0.02,
                                0.03, 0.04, 0.05}; // dont change etas :-)
    for (auto eta : etas) {
      auto [ops, e0] = triangular_12_complex(nup, eta);
      auto block = Spinhalf(nsites, nup);
      auto H = matrixC(ops, block, block);
      REQUIRE(H.is_hermitian(1e-8));

      arma::vec eigs;
      arma::eig_sym(eigs, H);

      // comment: reference data from Lanczos, only ~10 digits precise
      // Log("eigs(0): {}, e0: {}", eigs(0), e0);
      REQUIRE(std::abs(eigs(0) - e0) < 1e-8);
    }
  }

  {
    Log("spinhalf_matrix: Triangular J1J2Jchi N=12");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.toml";

    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    ops["J1"] = 1.00;
    ops["J2"] = 0.15;
    ops["Jchi"] = 0.09;

    int nsites = 12;
    int nup = 6;
    auto block = Spinhalf(nsites, nup);
    auto H = matrixC(ops, block, block);
    REQUIRE(H.is_hermitian(1e-8));

    arma::vec eigs;
    arma::eig_sym(eigs, H);
    double energy = -6.9456000700824329641;

    // Log("{:.18f} {:.18f}", eigs(0), energy);

    REQUIRE(isapprox(eigs(0), energy));
  }

  // Test S+/S-/Sz
  {
    Log.out("spinhalf_matrix: sp sm commutator test");

    for (int nsites = 2; nsites < 5; ++nsites) {

      auto block_raw = Spinhalf(nsites);
      for (int nup = 1; nup < nsites; ++nup) {

        auto block = Spinhalf(nsites, nup);
        auto blockp = Spinhalf(nsites, nup + 1);
        auto blockm = Spinhalf(nsites, nup - 1);

        for (int i = 0; i < nsites; ++i)
          for (int j = 0; j < nsites; ++j) {

            OpSum sp_i_m;
            sp_i_m += Op("S+", i);
            auto sp_i_m_mat = matrix(sp_i_m, blockm, block);
            auto sp_i_mat = matrix(sp_i_m, block_raw, block_raw);

            OpSum sm_j_m;
            sm_j_m += Op("S-", j);
            auto sm_j_m_mat = matrix(sm_j_m, block, blockm);
            auto sm_j_mat = matrix(sm_j_m, block_raw, block_raw);

            OpSum sp_i_p;
            sp_i_p += Op("S+", i);
            auto sp_i_p_mat = matrix(sp_i_p, block, blockp);

            OpSum sm_j_p;
            sm_j_p += Op("S-", j);
            auto sm_j_p_mat = matrix(sm_j_p, blockp, block);

            auto C1 = sp_i_m_mat * sm_j_m_mat;
            auto C2 = sm_j_p_mat * sp_i_p_mat;
            arma::mat comm = C1 - C2;
            auto C1r = sp_i_mat * sm_j_mat;
            auto C2r = sm_j_mat * sp_i_mat;
            arma::mat commr = C1r - C2r;

            if (i == j) {
              OpSum sz;
              sz += Op("Sz", i);
              auto sz_mat = matrix(sz, block, block);
              auto sz_matr = matrix(sz, block_raw, block_raw);
              REQUIRE(isapprox(comm, arma::mat(2.0 * sz_mat)));
              REQUIRE(isapprox(commr, arma::mat(2.0 * sz_matr)));
            } else {
              REQUIRE(isapprox(comm, arma::mat(comm.n_rows, comm.n_cols,
                                               arma::fill::zeros)));
              REQUIRE(isapprox(commr, arma::mat(commr.n_rows, commr.n_cols,
                                                arma::fill::zeros)));
            }
          }
      }
    }
  }

  test_onsite("Sz", "SzSz");

  for (int nsites = 2; nsites < 5; ++nsites) {
    auto b = Spinhalf(nsites);
    for (int s = 0; s < nsites; ++s) {

      // Exchange
      arma::mat m1 = 0.5 * (matrix(Op("S+", s), b) * matrix(Op("S-", s), b) +
                            matrix(Op("S-", s), b) * matrix(Op("S+", s), b));
      arma::mat m2 = matrix(Op("Exchange", {s, s}), b);
      REQUIRE(isapprox(m1, m2));

      // SdotS
      m1 = 0.5 * (matrix(Op("S+", s), b) * matrix(Op("S-", s), b) +
                  matrix(Op("S-", s), b) * matrix(Op("S+", s), b)) +
           matrix(Op("Sz", s), b) * matrix(Op("Sz", s), b);
      m2 = matrix(Op("SdotS", {s, s}), b);
      REQUIRE(isapprox(m1, m2));
    }
  }

} catch (xdiag::Error const &e) {
  xdiag::error_trace(e);
}
