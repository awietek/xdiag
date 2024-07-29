#include "../../catch.hpp"
#include <xdiag/basis/tj/basis_tj.hpp>
#include <xdiag/blocks/tj.hpp>

#include <xdiag/common.hpp>

#include "../../blocks/electron/testcases_electron.hpp"

template <typename bit_t>
void test_iterator_tj_basis_np(int64_t n_sites, int64_t n_up, int64_t n_dn) {
  using namespace xdiag;
  using basis_t = xdiag::basis::tj::BasisNp<bit_t>;
  basis_t basis(n_sites, n_up, n_dn);

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
    auto block = tJ(n_sites, n_up, n_dn);
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
void test_iterator_tj_basis_symmetric_np(int64_t n_sites, int64_t n_up,
                                         int64_t n_dn) {
  using namespace xdiag;

  auto [group, irreps] = testcases::electron::get_cyclic_group_irreps(n_sites);

  for (auto irrep : irreps) {
    auto basis =
        basis::tj::BasisSymmetricNp<bit_t>(n_sites, n_up, n_dn, group, irrep);
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
      auto block = tJ(n_sites, n_up, n_dn, group, irrep);
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
  for (int64_t n_sites = 0; n_sites <= 6; ++n_sites) {
    for (int64_t n_up = 0; n_up <= n_sites; ++n_up) {
      for (int64_t n_dn = 0; n_dn <= n_sites - n_up; ++n_dn) {
        test_iterator_tj_basis_np<uint32_t>(n_sites, n_up, n_dn);
        test_iterator_tj_basis_np<uint64_t>(n_sites, n_up, n_dn);
      }
    }
  }

  Log("Test tJ BasisSymmetricNp");
  for (int64_t n_sites = 0; n_sites <= 6; ++n_sites) {
    for (int64_t n_up = 0; n_up <= n_sites; ++n_up) {
      for (int64_t n_dn = 0; n_dn <= n_sites - n_up; ++n_dn) {
        test_iterator_tj_basis_symmetric_np<uint32_t>(n_sites, n_up, n_dn);
        test_iterator_tj_basis_symmetric_np<uint64_t>(n_sites, n_up, n_dn);
      }
    }
  }

} catch (xdiag::Error e) {
  error_trace(e);
}
