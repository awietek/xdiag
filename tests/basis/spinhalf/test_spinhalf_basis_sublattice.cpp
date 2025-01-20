#include "../../catch.hpp"

#include <iostream>

#include <xdiag/basis/spinhalf/basis_sublattice.hpp>
#include <xdiag/basis/spinhalf/basis_symmetric_no_sz.hpp>
#include <xdiag/basis/spinhalf/basis_symmetric_sz.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/common.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>

using namespace xdiag;

template <typename bit_t, class Basis1, class Basis2>
void compare_state(bit_t state, Basis1 const &basis1, Basis2 const &basis2) {
  auto ga1 = basis1.group_action();
  auto ga2 = basis2.group_action();
  int64_t idx1 = basis1.index(state);
  int64_t idx2 = basis2.index(state);
  if (idx1 != invalid_index) {
    REQUIRE(idx1 == idx2);
    bit_t rep1 = basis1.representative(state);
    bit_t rep2 = basis2.representative(state);
    REQUIRE(rep1 == rep2);
    bit_t norm1 = basis1.norm(idx1);
    bit_t norm2 = basis2.norm(idx2);
    REQUIRE(norm1 == norm2);
    auto [idx1a, sym1] = basis1.index_sym(state);
    auto [idx2a, sym2] = basis2.index_sym(state);
    REQUIRE(idx1a == idx1);
    REQUIRE(idx2a == idx2);
    REQUIRE(ga1.apply(sym1, state) == rep1);
    REQUIRE(ga2.apply(sym2, state) == rep2);
    auto [idx1b, syms1] = basis1.index_syms(state);
    auto [idx2b, syms2] = basis2.index_syms(state);
    REQUIRE(idx1b == idx1);
    REQUIRE(idx2b == idx2);
    REQUIRE(syms1.size() > 0);
    REQUIRE(syms2.size() > 0);
    for (auto s : syms1) {
      REQUIRE(ga1.apply(s, state) == rep1);
    }
    for (auto s : syms2) {
      REQUIRE(ga2.apply(s, state) == rep2);
    }
  }
}

template <typename bit_t, class Basis1, class Basis2>
void compare_indices_sz(Basis1 const &basis1, Basis2 const &basis2, int nup) {
  REQUIRE(basis1.size() == basis2.size());
  for (int64_t idx = 0; idx < basis1.size(); ++idx) {
    REQUIRE(basis1.state(idx) == basis2.state(idx));
  }

  int n_sites = basis1.n_sites();
  REQUIRE(n_sites == basis2.n_sites());
  for (auto state : combinatorics::Combinations<bit_t>(n_sites, nup)) {
    compare_state(state, basis1, basis2);
  }
}

template <typename bit_t, class Basis1, class Basis2>
void compare_indices_no_sz(Basis1 const &basis1, Basis2 const &basis2) {
  REQUIRE(basis1.size() == basis2.size());
  for (int64_t idx = 0; idx < basis1.size(); ++idx) {
    REQUIRE(basis1.state(idx) == basis2.state(idx));
  }

  int n_sites = basis1.n_sites();
  REQUIRE(n_sites == basis2.n_sites());
  for (auto state : combinatorics::Subsets<bit_t>(n_sites)) {
    compare_state(state, basis1, basis2);
  }
}

template <class bit_t> void test_spinhalf_basis_sublattice() {
  using basis_no_sz_t = basis::spinhalf::BasisSymmetricNoSz<bit_t>;
  using basis_sz_t = basis::spinhalf::BasisSymmetricSz<bit_t>;
  using basis_sl1_t = basis::spinhalf::BasisSublattice<bit_t, 1>;
  using basis_sl2_t = basis::spinhalf::BasisSublattice<bit_t, 2>;
  using basis_sl3_t = basis::spinhalf::BasisSublattice<bit_t, 3>;
  using basis_sl4_t = basis::spinhalf::BasisSublattice<bit_t, 4>;
  using basis_sl5_t = basis::spinhalf::BasisSublattice<bit_t, 5>;

  {
    Log("basis_spinhalf_sublattice: 1 sublattice");
    int n_sites = 8;
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.8.heisenberg.2sl.toml";
    auto fl = FileToml(lfile);
    std::vector<std::string> irrep_names = {
        "Gamma.D4.A1", "Gamma.D4.A2", "Gamma.D4.B1", "Gamma.D4.B2",
        "Gamma.D4.E",  "M.D4.A1",     "M.D4.A2",     "M.D4.B1",
        "M.D4.B2",     "M.D4.E",      "Sigma.D1.A",  "Sigma.D1.B",
        "X.D2.A1",     "X.D2.A2",     "X.D2.B1",     "X.D2.B2"};

    for (auto irrep_name : irrep_names) {
      auto irrep = read_representation(fl, irrep_name);
      auto idxng = basis_no_sz_t(irrep);
      auto idxng_sl = basis_sl1_t(irrep);

      compare_indices_no_sz<bit_t, basis_no_sz_t, basis_sl1_t>(idxng, idxng_sl);

      for (int nup = 0; nup <= n_sites; ++nup) {
        auto idxng = basis_sz_t(nup, irrep);
        auto idxng_sl = basis_sl1_t(nup, irrep);
        compare_indices_sz<bit_t, basis_sz_t, basis_sl1_t>(idxng, idxng_sl,
                                                           nup);
      }
    }
  }

  {
    Log("basis_spinhalf_sublattice: 2 sublattice");
    int n_sites = 8;
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.8.heisenberg.2sl.toml";
    auto fl = FileToml(lfile);
    std::vector<std::string> irrep_names = {
        "Gamma.D4.A1", "Gamma.D4.A2", "Gamma.D4.B1", "Gamma.D4.B2",
        "Gamma.D4.E",  "M.D4.A1",     "M.D4.A2",     "M.D4.B1",
        "M.D4.B2",     "M.D4.E",      "Sigma.D1.A",  "Sigma.D1.B",
        "X.D2.A1",     "X.D2.A2",     "X.D2.B1",     "X.D2.B2"};

    for (auto irrep_name : irrep_names) {
      auto irrep = read_representation(fl, irrep_name);
      auto idxng = basis_no_sz_t(irrep);
      auto idxng_sl = basis_sl2_t(irrep);
      compare_indices_no_sz<bit_t, basis_no_sz_t, basis_sl2_t>(idxng, idxng_sl);
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto idxng = basis_sz_t(nup, irrep);
        auto idxng_sl = basis_sl2_t(nup, irrep);
        compare_indices_sz<bit_t, basis_sz_t, basis_sl2_t>(idxng, idxng_sl,
                                                           nup);
      }
    }
  }

  {
    Log("basis_spinhalf_sublattice: 3 sublattice");
    int n_sites = 9;
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.9.heisenberg.3sl.toml";
    auto fl = FileToml(lfile);
    std::vector<std::string> irrep_names = {
        "Gamma.D2.A1", "Gamma.D2.A2", "Gamma.D2.B1",
        "Gamma.D2.B2", "Delta.C1.A",  "Sigma0.D1.A",
        "Sigma0.D1.B", "Sigma1.D1.A", "Sigma1.D1.B"};

    for (auto irrep_name : irrep_names) {
      auto irrep = read_representation(fl, irrep_name);
      auto idxng = basis_no_sz_t(irrep);
      auto idxng_sl = basis_sl3_t(irrep);
      compare_indices_no_sz<bit_t, basis_no_sz_t, basis_sl3_t>(idxng, idxng_sl);
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto idxng = basis_sz_t(nup, irrep);
        auto idxng_sl = basis_sl3_t(nup, irrep);
        compare_indices_sz<bit_t, basis_sz_t, basis_sl3_t>(idxng, idxng_sl,
                                                           nup);
      }
    }
  }

  {
    Log("basis_spinhalf_sublattice: 3 sublattice (triangular)");
    int n_sites = 9;
    std::string lfile = XDIAG_DIRECTORY
        "/misc/data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.toml";
    auto fl = FileToml(lfile);
    std::vector<std::string> irrep_names = {
        "Gamma.D6.A1", "Gamma.D6.A2", "Gamma.D6.B1", "Gamma.D6.B2",
        "Gamma.D6.E1", "Gamma.D6.E2", "K.D3.A1",     "K.D3.A2",
        "K.D3.E",      "Y.D1.A",      "Y.D1.B"};

    for (auto irrep_name : irrep_names) {
      auto irrep = read_representation(fl, irrep_name);
      auto idxng = basis_no_sz_t(irrep);
      auto idxng_sl = basis_sl3_t(irrep);
      compare_indices_no_sz<bit_t, basis_no_sz_t, basis_sl3_t>(idxng, idxng_sl);
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto idxng = basis_sz_t(nup, irrep);
        auto idxng_sl = basis_sl3_t(nup, irrep);
        compare_indices_sz<bit_t, basis_sz_t, basis_sl3_t>(idxng, idxng_sl,
                                                           nup);
      }
    }
  }

  {
    Log("basis_spinhalf_sublattice: 4 sublattice");
    int n_sites = 8;
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.8.heisenberg.4sl.toml";
    auto fl = FileToml(lfile);
    std::vector<std::string> irrep_names = {
        "Gamma.D4.A1", "Gamma.D4.A2", "Gamma.D4.B1", "Gamma.D4.B2",
        "Gamma.D4.E",  "M.D4.A1",     "M.D4.A2",     "M.D4.B1",
        "M.D4.B2",     "M.D4.E",      "Sigma.D1.A",  "Sigma.D1.B",
        "X.D2.A1",     "X.D2.A2",     "X.D2.B1",     "X.D2.B2"};

    for (auto irrep_name : irrep_names) {
      auto irrep = read_representation(fl, irrep_name);
      auto idxng = basis_no_sz_t(irrep);
      auto idxng_sl = basis_sl4_t(irrep);
      compare_indices_no_sz<bit_t, basis_no_sz_t, basis_sl4_t>(idxng, idxng_sl);
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto idxng = basis_sz_t(nup, irrep);
        auto idxng_sl = basis_sl4_t(nup, irrep);
        compare_indices_sz<bit_t, basis_sz_t, basis_sl4_t>(idxng, idxng_sl,
                                                           nup);
      }
    }
  }

  {
    Log("basis_spinhalf_sublattice: 5 sublattice");
    int n_sites = 10;
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.10.heisenberg.5sl.toml";
    auto fl = FileToml(lfile);
    std::vector<std::string> irrep_names = {
        "Gamma.C2.A", "Gamma.C2.B", "Delta0.C1.A", "Delta1.C1.A", "X.C2.A",
        "X.C2.B",     "Z0.C1.A",    "Z1.C1.A",     "Z2.C1.A",     "Z3.C1.A"};

    for (auto irrep_name : irrep_names) {
      auto irrep = read_representation(fl, irrep_name);
      auto idxng = basis_no_sz_t(irrep);
      auto idxng_sl = basis_sl5_t(irrep);
      compare_indices_no_sz<bit_t, basis_no_sz_t, basis_sl5_t>(idxng, idxng_sl);
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto idxng = basis_sz_t(nup, irrep);
        auto idxng_sl = basis_sl5_t(nup, irrep);
        compare_indices_sz<bit_t, basis_sz_t, basis_sl5_t>(idxng, idxng_sl,
                                                           nup);
      }
    }
  }
}

TEST_CASE("basis_spinhalf_sublattice", "[symmetries]") {
  Log("Test basis_spinhalf_sublattice");
  Log("uint32_t");
  test_spinhalf_basis_sublattice<uint32_t>();
  Log("uint64_t");
  test_spinhalf_basis_sublattice<uint64_t>();
  Log("Done");
}
