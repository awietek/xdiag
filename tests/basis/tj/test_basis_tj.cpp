// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"
#include <xdiag/basis/tj/basis_tj.hpp>
#include <xdiag/blocks/tj.hpp>

#include <xdiag/common.hpp>

#include "../../blocks/electron/testcases_electron.hpp"

template <typename bit_t>
void test_iterator_tj_basis_np(int64_t nsites, int64_t nup, int64_t ndn) {
  using namespace xdiag;
  using basis_t = xdiag::basis::tj::BasisNp<bit_t>;
  basis_t basis(nsites, nup, ndn);

  bit_t ups_prev = 0;
  bit_t dns_prev = 0;
  int64_t idx = 0;
  for (auto [ups, dns] : basis) {
    if (idx != 0) {
      if (ups == ups_prev) {
        REQUIRE(dns > dns_prev);
        dns_prev = dns;
      } else {
        REQUIRE(ups > ups_prev);
        ups_prev = ups;
        dns_prev = dns;
      }
    }
    ++idx;
  }
  REQUIRE(idx == basis.dim());

  {
    auto block = tJ(nsites, nup, ndn);
    int64_t idx = 0;
    for (auto pstate : block) {
      int64_t idx2 = block.index(pstate);
      // Log("{} {} {}", to_string(pstate), idx, idx2);
      REQUIRE(idx == idx2);
      ++idx;
    }
    REQUIRE(idx == block.dim());
  }
}

template <typename bit_t>
void test_iterator_tj_basis_symmetric_np(int64_t nsites, int64_t nup,
                                         int64_t ndn) {
  using namespace xdiag;

  auto irreps = testcases::electron::get_cyclic_group_irreps(nsites);

  for (auto irrep : irreps) {

    if (isreal(irrep)) {
      Vector characters = irrep.characters();
      auto basis = basis::tj::BasisSymmetricNp<bit_t>(
          nsites, nup, ndn, irrep.group(), characters.as<arma::vec>());
      bit_t ups_prev = 0;
      bit_t dns_prev = 0;
      int64_t idx = 0;
      for (auto [ups, dns] : basis) {
        if (idx != 0) {
          if (ups == ups_prev) {
            REQUIRE(dns > dns_prev);
            dns_prev = dns;
          } else {
            REQUIRE(ups > ups_prev);
            ups_prev = ups;
            dns_prev = dns;
          }
        }
        ++idx;
      }
      REQUIRE(idx == basis.dim());
    } else {
      Vector characters = irrep.characters();
      auto basis = basis::tj::BasisSymmetricNp<bit_t>(
          nsites, nup, ndn, irrep.group(), characters.as<arma::cx_vec>());
      bit_t ups_prev = 0;
      bit_t dns_prev = 0;
      int64_t idx = 0;
      for (auto [ups, dns] : basis) {
        if (idx != 0) {
          if (ups == ups_prev) {
            REQUIRE(dns > dns_prev);
            dns_prev = dns;
          } else {
            REQUIRE(ups > ups_prev);
            ups_prev = ups;
            dns_prev = dns;
          }
        }
        ++idx;
      }
      REQUIRE(idx == basis.dim());
    }

    {
      auto block = tJ(nsites, nup, ndn, irrep);
      int64_t idx = 0;
      for (auto pstate : block) {
        int64_t idx2 = block.index(pstate);
        // Log("{} {} {}", to_string(pstate), idx, idx2);
        REQUIRE(idx == idx2);
        ++idx;
      }
      REQUIRE(idx == block.dim());
    }
  }
}

TEST_CASE("tj_basis", "[basis]") try {
  using namespace xdiag;

  Log("Test tJ BasisNp");
  for (int64_t nsites = 1; nsites <= 6; ++nsites) {
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; ndn <= nsites - nup; ++ndn) {
        test_iterator_tj_basis_np<uint32_t>(nsites, nup, ndn);
        test_iterator_tj_basis_np<uint64_t>(nsites, nup, ndn);
      }
    }
  }

  Log("Test tJ BasisSymmetricNp");
  for (int64_t nsites = 1; nsites <= 6; ++nsites) {
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; ndn <= nsites - nup; ++ndn) {
        test_iterator_tj_basis_symmetric_np<uint32_t>(nsites, nup, ndn);
        test_iterator_tj_basis_symmetric_np<uint64_t>(nsites, nup, ndn);
      }
    }
  }

} catch (xdiag::Error e) {
  error_trace(e);
}
