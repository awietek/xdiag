#pragma once

#include <xdiag/common.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/bondlist.hpp>
#include <xdiag/operators/couplings.hpp>

#include <lila/all.hpp>

namespace xdiag::electron {

template <class bit_t, class coeff_t, class Filler>
void do_U(Couplings const &couplings, Electron<bit_t> const &block,
          lila::Vector<coeff_t> const &vec_in, lila::Vector<coeff_t> &vec_out) {


  // if (couplings.defined("U") && !lila::close(couplings["U"], (complex)0.)) {

  //   int nsites = block.nsites();
  //   int ndn = block.ndn();

  //   double U = lila::real(couplings["U"]);

  //   int64_t idx = 0;
  //   for (bit_t ups : my_ups_) // loop over upspins of process
  //   {
  //     for (bit_t dns :
  //          Combinations<bit_t>(nsites, ndn)) // loop over all downspins
  //     {
  //       double coeff = U * (double)popcnt(ups & dns);
  //       vec_out(idx) += coeff * vec_in(idx);
  //       ++idx;
  //     }
  //   }
  // }
}

} // namespace xdiag::electron
