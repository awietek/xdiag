#include "../catch.hpp"

#include <iostream>

#include <xdiag/common.hpp>
#include <xdiag/symmetries/representation.hpp>
#include "../blocks/electron/testcases_electron.hpp"
#include <xdiag/utils/xdiag_show.hpp>

using namespace xdiag;
using namespace arma;
// using namespace std::complex_literals;

TEST_CASE("representation", "[symmetries]") try {
//   Log("Testing real representation construction");
  for (int64_t n_sites = 3; n_sites < 8; ++n_sites){
    auto irreps = testcases::electron::get_cyclic_group_irreps(n_sites);
    for (auto irrep:irreps) {
        auto chars_hc = arma::conj(irrep.characters().as<arma::cx_vec>());
        auto irrep_hc = Representation(irrep.group(), arma::cx_vec(chars_hc));
        // XDIAG_SHOW(irrep);
        // XDIAG_SHOW(irrep_hc);
        auto irrep_id = irrep * irrep_hc;
        // XDIAG_SHOW(irrep_id);
        REQUIRE(irrep_id.isreal());        
    }
  }
} catch (xdiag::Error e) {
    xdiag::error_trace(e);
}
