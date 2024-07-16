#include "../../catch.hpp"

#include "testcases_spinhalf.hpp"
#include <iostream>
#include <xdiag/blocks/spinhalf/spinhalf_matrix.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/utils/close.hpp>
#include <xdiag/utils/print_macro.hpp>

using namespace xdiag;

TEST_CASE("spinhalf_matrix", "[spinhalf]") {
  using namespace xdiag::testcases::spinhalf;

  {
    Log.out("spinhalf_matrix: Heisenberg chain test, J=1.0, N=2,..,6");
    for (int n_sites = 2; n_sites <= 6; ++n_sites) {
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto [ops, exact_eigs] = HBchain_fullspectrum_nup(n_sites, nup);
        auto block = Spinhalf(n_sites, nup);
        auto H = matrix(ops, block, block);
        REQUIRE(H.is_hermitian(1e-7));
        arma::vec eigs;
        arma::eig_sym(eigs, H);
        REQUIRE(close(eigs, exact_eigs));
      }
    }
  }

  {
    Log.out("spinhalf_matrix: Heisenberg all-to-all tJ comparison");
    for (int n_sites = 2; n_sites <= 6; ++n_sites)
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto ops = HB_alltoall(n_sites);
        auto block = Spinhalf(n_sites, nup);
        auto block_tJ = tJ(n_sites, nup, n_sites - nup);
        auto H = matrix(ops, block, block);
        auto H_tJ = matrix(ops, block_tJ, block_tJ);
        REQUIRE(H.is_hermitian());
        REQUIRE(H_tJ.is_hermitian());

        arma::vec eigs;
        arma::eig_sym(eigs, H);

        arma::vec eigs_tJ;
        arma::eig_sym(eigs_tJ, H_tJ);
        REQUIRE(close(eigs, eigs_tJ));
      }
  }

  {
    Log.out("spinhalf_matrix: Heisenberg all-to-all Sz <-> NoSz comparison");
    for (int n_sites = 2; n_sites <= 6; ++n_sites) {
      auto ops = HB_alltoall(n_sites);
      auto block_no_sz = Spinhalf(n_sites);
      auto H_no_sz = matrix(ops, block_no_sz, block_no_sz);
      REQUIRE(H_no_sz.is_hermitian(1e-8));
      arma::vec eigs_no_sz;
      arma::eig_sym(eigs_no_sz, H_no_sz);

      std::vector<double> eigs_sz_all;

      for (int nup = 0; nup <= n_sites; ++nup) {
        auto block_sz = Spinhalf(n_sites, nup);
        auto H_sz = matrix(ops, block_sz, block_sz);
        REQUIRE(H_sz.is_hermitian(1e-7));
        arma::vec eigs_sz;
        arma::eig_sym(eigs_sz, H_sz);

        for (auto eig : eigs_sz)
          eigs_sz_all.push_back(eig);
      }
      std::sort(eigs_sz_all.begin(), eigs_sz_all.end());

      REQUIRE(close(eigs_no_sz, arma::vec(eigs_sz_all)));
    }
  }

  {
    Log.out("spinhalf_matrix: triangular N=12 complex exchange");
    int n_sites = 12;
    int nup = 6;
    std::vector<double> etas = {0.00, 0.01, 0.02,
                                0.03, 0.04, 0.05}; // dont change etas :-)
    for (auto eta : etas) {
      auto [ops, e0] = triangular_12_complex(nup, eta);
      auto block = Spinhalf(n_sites, nup);
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
                        "triangular.12.j1j2jch.sublattices.fsl.lat";

    auto ops = read_opsum(lfile);
    ops["J1"] = 1.00;
    ops["J2"] = 0.15;
    ops["Jchi"] = 0.09;

    int n_sites = 12;
    int n_up = 6;
    auto block = Spinhalf(n_sites, n_up);
    auto H = matrixC(ops, block, block);
    REQUIRE(H.is_hermitian(1e-8));

    arma::vec eigs;
    arma::eig_sym(eigs, H);
    double energy = -6.9456000700824329641;

    // Log("{:.18f} {:.18f}", eigs(0), energy);

    REQUIRE(close(eigs(0), energy));
  }

  // Test S+/S-/Sz
  {
    Log.out("spinhalf_matrix: sp sm commutator test");

    for (int n_sites = 2; n_sites < 5; ++n_sites) {

      auto block_raw = Spinhalf(n_sites);
      for (int nup = 1; nup < n_sites; ++nup) {

        auto block = Spinhalf(n_sites, nup);
        auto blockp = Spinhalf(n_sites, nup + 1);
        auto blockm = Spinhalf(n_sites, nup - 1);

        for (int i = 0; i < n_sites; ++i)
          for (int j = 0; j < n_sites; ++j) {

            OpSum sp_i_m;
            sp_i_m += Op("S+", 1.0, i);
            auto sp_i_m_mat = matrix(sp_i_m, blockm, block);
            auto sp_i_mat = matrix(sp_i_m, block_raw, block_raw);

            OpSum sm_j_m;
            sm_j_m += Op("S-", 1.0, j);
            auto sm_j_m_mat = matrix(sm_j_m, block, blockm);
            auto sm_j_mat = matrix(sm_j_m, block_raw, block_raw);

            OpSum sp_i_p;
            sp_i_p += Op("S+", 1.0, i);
            auto sp_i_p_mat = matrix(sp_i_p, block, blockp);

            OpSum sm_j_p;
            sm_j_p += Op("S-", 1.0, j);
            auto sm_j_p_mat = matrix(sm_j_p, blockp, block);

            auto C1 = sp_i_m_mat * sm_j_m_mat;
            auto C2 = sm_j_p_mat * sp_i_p_mat;
            arma::mat comm = C1 - C2;
            auto C1r = sp_i_mat * sm_j_mat;
            auto C2r = sm_j_mat * sp_i_mat;
            arma::mat commr = C1r - C2r;

            if (i == j) {
              OpSum sz;
              sz += Op("SZ", 1.0, i);
              auto sz_mat = matrix(sz, block, block);
              auto sz_matr = matrix(sz, block_raw, block_raw);
              REQUIRE(close(comm, arma::mat(2.0 * sz_mat)));
              REQUIRE(close(commr, arma::mat(2.0 * sz_matr)));
            } else {
              REQUIRE(close(comm, arma::mat(comm.n_rows, comm.n_cols,
                                            arma::fill::zeros)));
              REQUIRE(close(commr, arma::mat(commr.n_rows, commr.n_cols,
                                             arma::fill::zeros)));
            }
          }
      }
    }
  }
}
