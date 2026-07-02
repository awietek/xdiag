// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <iostream>

#include <xdiag/basis/basis_sublattice.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/config.hpp>

using namespace xdiag;

// Check internal consistency of a BasisSublattice state lookup
template <typename bit_t, int n_sublat>
void check_state(bit_t state,
                 basis::BasisSublattice<bit_t, n_sublat> const &basis) {
  auto const &action = basis.action();
  int64_t idx = basis.index(state);
  if (idx == -1)
    return;

  bit_t rep = basis.representative(state);
  REQUIRE(basis[idx] == rep);

  auto [idx_sym, sym] = basis.index_sym(state);
  REQUIRE(idx_sym == idx);
  REQUIRE(action.apply(sym, state) == rep);

  auto [idx_syms, syms] = basis.index_syms(state);
  REQUIRE(idx_syms == idx);
  REQUIRE(syms.size() > 0);
  for (auto s : syms) {
    REQUIRE(action.apply(s, state) == rep);
  }
}

// Compare BasisSublattice against a BasisSymmetric reference
template <typename bit_t, int n_sublat, class BasisRef>
void compare_with_reference(
    basis::BasisSublattice<bit_t, n_sublat> const &basis_sl,
    BasisRef const &basis_ref) {
  REQUIRE(basis_sl.nsites() == basis_ref.nsites());
  REQUIRE(basis_sl.size() == basis_ref.size());
  for (int64_t idx = 0; idx < basis_sl.size(); ++idx) {
    REQUIRE(basis_sl[idx] == basis_ref[idx]);
    REQUIRE(basis_sl.norm(idx) == Approx(basis_ref.norm(idx)));
  }
}

template <typename bit_t, int n_sublat>
void test_sublattice_no_sz(Representation const &irrep) {
  using namespace basis;
  using namespace combinatorics;

  auto basis_sl =
      BasisSublattice<bit_t, n_sublat>(irrep.group(), irrep.characters());
  auto basis_ref = BasisSymmetric<Subsets<bit_t>>(
      Subsets<bit_t>(irrep.group().nsites()), irrep.group(),
      irrep.characters());

  compare_with_reference(basis_sl, basis_ref);

  int64_t nsites = basis_sl.nsites();
  for (auto state : Subsets<bit_t>(nsites)) {
    check_state(state, basis_sl);
  }
}

template <typename bit_t, int n_sublat>
void test_sublattice_sz(int64_t nup, Representation const &irrep) {
  using namespace basis;
  using namespace combinatorics;

  auto basis_sl = BasisSublattice<bit_t, n_sublat>(nup, irrep.group(),
                                                   irrep.characters());
  auto basis_ref = BasisSymmetric<Combinations<bit_t>>(
      Combinations<bit_t>(irrep.group().nsites(), nup), irrep.group(),
      irrep.characters());

  compare_with_reference(basis_sl, basis_ref);

  int64_t nsites = basis_sl.nsites();
  for (auto state : Combinations<bit_t>(nsites, nup)) {
    check_state(state, basis_sl);
  }
}

template <typename bit_t, int n_sublat>
void test_sublattice(int64_t nsites, std::string const &lfile,
                     std::vector<std::string> const &irrep_names) {
  auto fl = FileToml(lfile);
  for (auto const &irrep_name : irrep_names) {
    auto irrep = read_representation(fl, irrep_name);
    test_sublattice_no_sz<bit_t, n_sublat>(irrep);
    for (int nup = 0; nup <= nsites; ++nup) {
      test_sublattice_sz<bit_t, n_sublat>(nup, irrep);
    }
  }
}

template <typename bit_t> void test_basis_sublattice() {
  {
    Log("basis_sublattice: 1 sublattice");
    int nsites = 8;
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.8.heisenberg.2sl.toml";
    std::vector<std::string> irrep_names = {
        "Gamma.D4.A1", "Gamma.D4.A2", "Gamma.D4.B1", "Gamma.D4.B2",
        "Gamma.D4.E",  "M.D4.A1",     "M.D4.A2",     "M.D4.B1",
        "M.D4.B2",     "M.D4.E",      "Sigma.D1.A",  "Sigma.D1.B",
        "X.D2.A1",     "X.D2.A2",     "X.D2.B1",     "X.D2.B2"};
    test_sublattice<bit_t, 1>(nsites, lfile, irrep_names);
  }

  {
    Log("basis_sublattice: 2 sublattice");
    int nsites = 8;
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.8.heisenberg.2sl.toml";
    std::vector<std::string> irrep_names = {
        "Gamma.D4.A1", "Gamma.D4.A2", "Gamma.D4.B1", "Gamma.D4.B2",
        "Gamma.D4.E",  "M.D4.A1",     "M.D4.A2",     "M.D4.B1",
        "M.D4.B2",     "M.D4.E",      "Sigma.D1.A",  "Sigma.D1.B",
        "X.D2.A1",     "X.D2.A2",     "X.D2.B1",     "X.D2.B2"};
    test_sublattice<bit_t, 2>(nsites, lfile, irrep_names);
  }

  {
    Log("basis_sublattice: 3 sublattice (square)");
    int nsites = 9;
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.9.heisenberg.3sl.toml";
    std::vector<std::string> irrep_names = {
        "Gamma.D2.A1", "Gamma.D2.A2", "Gamma.D2.B1",
        "Gamma.D2.B2", "Delta.C1.A",  "Sigma0.D1.A",
        "Sigma0.D1.B", "Sigma1.D1.A", "Sigma1.D1.B"};
    test_sublattice<bit_t, 3>(nsites, lfile, irrep_names);
  }

  {
    Log("basis_sublattice: 3 sublattice (triangular)");
    int nsites = 9;
    std::string lfile = XDIAG_DIRECTORY
        "/misc/data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.toml";
    std::vector<std::string> irrep_names = {
        "Gamma.D6.A1", "Gamma.D6.A2", "Gamma.D6.B1", "Gamma.D6.B2",
        "Gamma.D6.E1", "Gamma.D6.E2", "K.D3.A1",     "K.D3.A2",
        "K.D3.E",      "Y.D1.A",      "Y.D1.B"};
    test_sublattice<bit_t, 3>(nsites, lfile, irrep_names);
  }

  {
    Log("basis_sublattice: 4 sublattice");
    int nsites = 8;
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.8.heisenberg.4sl.toml";
    std::vector<std::string> irrep_names = {
        "Gamma.D4.A1", "Gamma.D4.A2", "Gamma.D4.B1", "Gamma.D4.B2",
        "Gamma.D4.E",  "M.D4.A1",     "M.D4.A2",     "M.D4.B1",
        "M.D4.B2",     "M.D4.E",      "Sigma.D1.A",  "Sigma.D1.B",
        "X.D2.A1",     "X.D2.A2",     "X.D2.B1",     "X.D2.B2"};
    test_sublattice<bit_t, 4>(nsites, lfile, irrep_names);
  }

  {
    Log("basis_sublattice: 5 sublattice");
    int nsites = 10;
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.10.heisenberg.5sl.toml";
    std::vector<std::string> irrep_names = {
        "Gamma.C2.A", "Gamma.C2.B", "Delta0.C1.A", "Delta1.C1.A", "X.C2.A",
        "X.C2.B",     "Z0.C1.A",    "Z1.C1.A",     "Z2.C1.A",     "Z3.C1.A"};
    test_sublattice<bit_t, 5>(nsites, lfile, irrep_names);
  }
}

TEST_CASE("basis_sublattice", "[basis]") try {
  Log("Test basis_sublattice");
  Log("uint32_t");
  test_basis_sublattice<uint32_t>();
  Log("uint64_t");
  test_basis_sublattice<uint64_t>();
  Log("Done");
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}
