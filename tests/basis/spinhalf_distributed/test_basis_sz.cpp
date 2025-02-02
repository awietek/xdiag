#include "../../catch.hpp"

#include <mpi.h>

#include <xdiag/basis/spinhalf_distributed/basis_sz.hpp>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/parallel/mpi/allreduce.hpp>
#include <xdiag/utils/logger.hpp>

template <typename bit_t> void test_spinhalf_distributed_basis_sz() {
  using namespace xdiag;

  for (int nsites = 2; nsites <= 6; ++nsites) {
    for (int nup = 0; nup <= nsites; ++nup) {
      // int nsites=16;
      // int nup = 7;
      auto basis = basis::spinhalf_distributed::BasisSz<bit_t>(nsites, nup);

      // pre / post order
      bit_t spins_before = 0;
      int64_t idx = 0;
      // Log("N: {} nup: {}", nsites, nup);
      for (bit_t prefix : basis.prefixes()) {
        bit_t prefix_shifted = prefix << basis.n_postfix_bits();
        for (bit_t postfix : basis.postfix_states(prefix)) {
          bit_t spins = prefix_shifted | postfix;
          REQUIRE(bits::popcnt(spins) == nup);
          // std::cout << BSTR(spins) << "\n";
          if (idx > 0) {
            REQUIRE(spins > spins_before);
          }
          spins_before = spins;
          ++idx;
        }
      }
      REQUIRE(idx == basis.size());

      // post / pre order
      spins_before = 0;
      idx = 0;
      for (bit_t postfix : basis.postfixes()) {
        bit_t postfix_shifted = postfix << basis.n_prefix_bits();
        for (bit_t prefix : basis.prefix_states(postfix)) {
          bit_t spins = postfix_shifted | prefix;
          REQUIRE(bits::popcnt(spins) == nup);
          if (idx > 0) {
            REQUIRE(spins > spins_before);
          }
          spins_before = spins;
          ++idx;
        }
      }
      REQUIRE(idx == basis.size_transpose());
    }
  }
}

TEST_CASE("spinhalf_distributed_basis_sz", "[spinhalf_distributed]") {
  using namespace xdiag;

  Log("spinhalf_distributed_basis_sz test (uint32_t)");
  test_spinhalf_distributed_basis_sz<uint32_t>();
  Log("spinhalf_distributed_basis_sz test (uint64_t)");
  test_spinhalf_distributed_basis_sz<uint64_t>();
}
