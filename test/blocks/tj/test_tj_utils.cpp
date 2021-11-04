#include "../../catch.hpp"

#include <iostream>

#include "testcases_tj.h"
#include <hydra/all.h>
#include <lila/all.h>

using namespace hydra;

TEST_CASE("tj_utils", "[tj]") {
  using namespace hydra::combinatorics;
  using namespace hydra::utils;
  using namespace hydra::bitops;
  lila::Log.out("Testing tj_utils");

  for (int N = 0; N <= 6; ++N)
    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {
        // lila::Log.out("N: {},  nup: {}, ndn: {}", N, nup, ndn);

        for (auto up : Combinations(N, nup))
          for (auto dnc : Combinations(N - nup, ndn)) {
            auto [upc, dn] = up_dnc_to_upc_dn(up, dnc);
            auto [up2, dnc2] = upc_dn_to_up_dnc(upc, dn);
	    REQUIRE(up == up2);
	    REQUIRE(dnc == dnc2);
            // lila::Log.out("up : {} dnc: {}", bits_to_string(up, N), bits_to_string(dnc,N));
            // lila::Log.out("upc: {} dn : {}", bits_to_string(upc,N), bits_to_string(dn, N));
          }

        for (auto dn : Combinations(N, ndn))
          for (auto upc : Combinations(N - ndn, nup)) {
            auto [up, dnc] = upc_dn_to_up_dnc(upc, dn);
            auto [upc2, dn2] = up_dnc_to_upc_dn(up, dnc);
	    REQUIRE(upc == upc2);
	    REQUIRE(dn == dn2);
            // lila::Log.out("up : {} dnc: {}", bits_to_string(up, N), bits_to_string(dnc,N));
            // lila::Log.out("upc: {} dn : {}", bits_to_string(upc,N), bits_to_string(dn, N));
          }

      }
}
