#include "../../catch.hpp"
#include <xdiag/basis/electron/basis_electron.hpp>
#include <xdiag/blocks/electron.hpp>

#include <xdiag/common.hpp>
#include <xdiag/utils/logger.hpp>

#include "../../blocks/electron/testcases_electron.hpp"

template <typename bit_t>
void test_iterator_electron_basis_no_np(int64_t nsites) {
  using namespace xdiag;

  using basis_t = xdiag::basis::electron::BasisNoNp<bit_t>;
  basis_t basis(nsites);

  bit_t ups_prev = 0;
  bit_t dns_prev = 0;
  int64_t idx = 0;
  for (auto [ups, dns] : basis) {
    // Log("{} {} {} {} {}", idx, BSTR(ups), BSTR(dns), BSTR(ups_prev),
    // BSTR(dns_prev));
    if (idx == 0) {
      REQUIRE(ups == ups_prev);
      REQUIRE(dns == dns_prev);
    } else {
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
    auto block = Electron(nsites);
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
void test_iterator_electron_basis_symmetric_no_np(int64_t nsites) {
  using namespace xdiag;

  auto irreps = xdiag::testcases::electron::get_cyclic_group_irreps(nsites);

  for (auto irrep : irreps) {
    auto basis =
        xdiag::basis::electron::BasisSymmetricNoNp<bit_t>(nsites, irrep);

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
      auto block = Electron(nsites, irrep);
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

template <typename bit_t>
void test_iterator_electron_basis_np(int64_t nsites, int64_t nup,
                                     int64_t ndn) {
  using namespace xdiag;
  using basis_t = xdiag::basis::electron::BasisNp<bit_t>;
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
    auto block = Electron(nsites, nup, ndn);
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
void test_iterator_electron_basis_symmetric_np(int64_t nsites, int64_t nup,
                                               int64_t ndn) {
  using namespace xdiag;

  auto irreps = testcases::electron::get_cyclic_group_irreps(nsites);

  for (auto irrep : irreps) {
    auto basis =
        basis::electron::BasisSymmetricNp<bit_t>(nsites, nup, ndn, irrep);
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
      auto block = Electron(nsites, nup, ndn, irrep);
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

TEST_CASE("electron_basis", "[basis]") {
  using namespace xdiag;
  Log("Test Electron BasisNoNp");
  for (int64_t nsites = 1; nsites < 6; ++nsites) {
    test_iterator_electron_basis_no_np<uint32_t>(nsites);
    test_iterator_electron_basis_no_np<uint64_t>(nsites);
  }

  Log("Test Electron BasisNp");
  for (int64_t nsites = 1; nsites < 6; ++nsites) {
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
        test_iterator_electron_basis_np<uint32_t>(nsites, nup, ndn);
        test_iterator_electron_basis_np<uint64_t>(nsites, nup, ndn);
      }
    }
  }

  Log("Test Electron BasisSymmetricNoNp");
  for (int64_t nsites = 1; nsites < 6; ++nsites) {
    test_iterator_electron_basis_symmetric_no_np<uint32_t>(nsites);
    test_iterator_electron_basis_symmetric_no_np<uint64_t>(nsites);
  }

  Log("Test Electron BasisSymmetricNp");
  for (int64_t nsites = 1; nsites < 6; ++nsites) {
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
        test_iterator_electron_basis_symmetric_np<uint32_t>(nsites, nup,
                                                            ndn);
        test_iterator_electron_basis_symmetric_np<uint64_t>(nsites, nup,
                                                            ndn);
      }
    }
  }
}
