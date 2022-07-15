#include "../../catch.hpp"

#include "testcases_spinhalf.h"
#include <hydra/all.h>
#include <iostream>

using namespace hydra;

TEST_CASE("spinhalf_matrix", "[models][spinhalf]") {
  using namespace hydra::testcases::spinhalf;

  {
    Log.out("spinhalf_matrix: Heisenberg chain test, J=1.0, N=2,..,6");
    for (int n_sites = 2; n_sites <= 6; ++n_sites)
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto [bonds, couplings, exact_eigs] =
            HBchain_fullspectrum_nup(n_sites, nup);
        auto block = Spinhalf<uint32_t>(n_sites, nup);
        auto H = MatrixReal(bonds, couplings, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));
        auto eigs = lila::EigenvaluesSym(H);
        REQUIRE(lila::close(eigs, exact_eigs));
      }
  }

  {
    Log.out("spinhalf_matrix: Heisenberg all-to-all tJ comparison");
    for (int n_sites = 2; n_sites <= 6; ++n_sites)
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto [bonds, couplings] = HB_alltoall(n_sites);

        auto block = Spinhalf<uint32_t>(n_sites, nup);
        auto block_tJ = tJ<uint32_t>(n_sites, nup, n_sites - nup);
        auto H = MatrixReal(bonds, couplings, block, block);
        auto H_tJ = MatrixReal(bonds, couplings, block_tJ, block_tJ);
        REQUIRE(lila::close(H, lila::Herm(H)));
        REQUIRE(lila::close(H_tJ, lila::Herm(H_tJ)));

        auto eigs = lila::EigenvaluesSym(H);
        auto eigs_tJ = lila::EigenvaluesSym(H_tJ);
        REQUIRE(lila::close(eigs, eigs_tJ));
      }
  }

  {
    Log.out(
        "spinhalf_matrix: Heisenberg all-to-all Sz <-> NoSz comparison");
    for (int n_sites = 2; n_sites <= 6; ++n_sites) {
      auto [bonds, couplings] = HB_alltoall(n_sites);

      auto block_no_sz = Spinhalf(n_sites);
      auto H_no_sz = MatrixReal(bonds, couplings, block_no_sz, block_no_sz);
      REQUIRE(lila::close(H_no_sz, lila::Herm(H_no_sz)));
      auto eigs_no_sz = lila::EigenvaluesSym(H_no_sz);

      lila::Vector<double> eigs_sz_all;

      for (int nup = 0; nup <= n_sites; ++nup) {
        auto block_sz = Spinhalf(n_sites, nup);
        auto H_sz = MatrixReal(bonds, couplings, block_sz, block_sz);
        REQUIRE(lila::close(H_sz, lila::Herm(H_sz)));
        auto eigs_sz = lila::EigenvaluesSym(H_sz);

        for (auto eig : eigs_sz)
          eigs_sz_all.push_back(eig);
      }
      std::sort(eigs_sz_all.begin(), eigs_sz_all.end());

      REQUIRE(lila::close(eigs_no_sz, eigs_sz_all));
    }
  }

  {
    Log.out("spinhalf_matrix: triangular N=12 complex exchange");
    int n_sites = 12;
    int nup = 6;
    std::vector<double> etas = {0.00, 0.01, 0.02,
                                0.03, 0.04, 0.05}; // dont change etas :-)
    for (auto eta : etas) {
      auto [bonds, couplings, e0] = triangular_12_complex(nup, eta);
      auto block = Spinhalf<uint32_t>(n_sites, nup);
      auto H = MatrixCplx(bonds, couplings, block, block);
      REQUIRE(lila::close(H, lila::Herm(H)));

      auto eigs = lila::EigenvaluesSym(H);

      // comment: reference data from Lanczos, only ~10 digits precise
      // Log("eigs(0): {}, e0: {}", eigs(0), e0);
      REQUIRE(std::abs(eigs(0) - e0) < 1e-8);
    }
  }

  {
    Log("spinhalf_matrix: Triangular J1J2Jchi N=12");
    std::string lfile = "data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.lat";

    auto bondlist = read_bondlist(lfile);
    Couplings couplings;
    couplings["J1"] = 1.00;
    couplings["J2"] = 0.15;
    couplings["Jchi"] = 0.09;

    int n_sites = 12;
    int n_up = 6;
    auto block = Spinhalf<uint16_t>(n_sites, n_up);
    auto H = MatrixCplx(bondlist, couplings, block, block);
    REQUIRE(lila::close(H, lila::Herm(H)));
    
    auto eigs = lila::EigenvaluesSym(H);
    double energy = -6.9456000700824329641;
    
    // Log("{:.18f} {:.18f}", eigs(0), energy);

    REQUIRE(lila::close(eigs(0), energy));

  }

  // Test S+/S-/Sz
  {
    Log.out("spinhalf_matrix: sp sm commutator test");

    for (int n_sites = 2; n_sites < 5; ++n_sites) {

      Couplings cpls;
      cpls["H"] = 1.0;

      auto block_raw = Spinhalf(n_sites);
      for (int nup = 1; nup < n_sites; ++nup) {

	auto block = Spinhalf(n_sites, nup);
	auto blockp = Spinhalf(n_sites, nup + 1);
	auto blockm = Spinhalf(n_sites, nup - 1);

        for (int i = 0; i < n_sites; ++i)
          for (int j = 0; j < n_sites; ++j) {


            BondList sp_i_m;
            sp_i_m << Bond("S+", "H", i);
            auto sp_i_m_mat = MatrixReal(sp_i_m, cpls, blockm, block);
            auto sp_i_mat = MatrixReal(sp_i_m, cpls, block_raw, block_raw);
	   
            BondList sm_j_m;
            sm_j_m << Bond("S-", "H", j);
            auto sm_j_m_mat = MatrixReal(sm_j_m, cpls, block, blockm);
            auto sm_j_mat = MatrixReal(sm_j_m, cpls, block_raw, block_raw);

            BondList sp_i_p;
            sp_i_p << Bond("S+", "H", i);
            auto sp_i_p_mat = MatrixReal(sp_i_p, cpls, block, blockp);

            BondList sm_j_p;
            sm_j_p << Bond("S-", "H", j);
            auto sm_j_p_mat = MatrixReal(sm_j_p, cpls, blockp, block);

            auto C1 = lila::Mult(sp_i_m_mat, sm_j_m_mat);
            auto C2 = lila::Mult(sm_j_p_mat, sp_i_p_mat);
            auto comm = C1 - C2;
            auto C1r = lila::Mult(sp_i_mat, sm_j_mat);
            auto C2r = lila::Mult(sm_j_mat, sp_i_mat);
            auto commr = C1r - C2r;

            if (i == j) {
              BondList sz;
              sz << Bond("SZ", "H", i);
              auto sz_mat = MatrixReal(sz, cpls, block, block);
              auto sz_matr = MatrixReal(sz, cpls, block_raw, block_raw);
              REQUIRE(lila::close(comm, 2.0 * sz_mat));
              REQUIRE(lila::close(commr, 2.0 * sz_matr));
            } else {
              REQUIRE(lila::close(comm, 0.0));
              REQUIRE(lila::close(commr, 0.0));
            }
          }
      }
    }
  }
}
