#pragma once

#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/blocks/spinhalf_mpi/spinhalf_mpi.h>
#include <hydra/mpi/logger_mpi.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>

namespace hydra::terms::spinhalf_mpi {

template <class bit_t, class coeff_t>
void do_ising_mpi(BondList const &bonds, Couplings const &couplings,
                  SpinhalfMPI<bit_t> const &block,
                  lila::Vector<coeff_t> const &vec_in,
                  lila::Vector<coeff_t> &vec_out) {

  auto ising = bonds.bonds_of_type("HEISENBERG") +
               bonds.bonds_of_type("ISING") + bonds.bonds_of_type("HB");

  for (auto bond : ising) {

    if (bond.size() != 2)
      LogMPI.err("Error computing SpinhalfMPI Ising: "
                 "bond must have exactly two sites defined");

    std::string coupling = bond.coupling();
    if (couplings.defined(coupling) &&
        !lila::close(couplings[coupling], (complex)0.)) {

      double J = lila::real(couplings[coupling]);

      // Set values for same/diff (tJ block definition)
      double val_same = J / 4.;
      double val_diff = -J / 4.;

      int s1 = std::min(bond.site(0), bond.site(1));
      int s2 = std::max(bond.site(0), bond.site(1));
      if (s1 == s2)
        LogMPI.err("Error computing SpinhalfMPI Ising: "
                   "operator acting on twice the same site");
      bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
      int n_postfix_bits = block.n_postfix_bits_;

      idx_t idx = 0;
      for (auto prefix : block.prefixes_) {
        int n_up_prefix = bitops::popcnt(prefix);
        int n_up_postfix = block.n_up() - n_up_prefix;
        bit_t prefix_shifted = (prefix << n_postfix_bits);

        auto const &postfixes = block.postfix_states_[n_up_postfix];
	// // Simple variant
        // for (auto postfix : postfixes) {
        //   bit_t state = prefix_shifted | postfix;
        //   double val = (bitops::popcnt(state & mask) == 1) ? val_diff :
        //   val_same; vec_out(idx) += val * vec_in(idx);
        //   ++idx;
        // }


	// Fast variant

        // Both sites are on prefixes
        if ((s1 >= n_postfix_bits) && (s2 >= n_postfix_bits)) {
          double val =
              (bitops::popcnt(prefix_shifted & mask) & 1) ? val_diff : val_same;
          idx_t end = idx + postfixes.size();
          for (; idx < end; ++idx) {
	    vec_out(idx) += val * vec_in(idx);
          }
        }

        // Both sites are on postfixes
        else if ((s1 < n_postfix_bits) && (s2 < n_postfix_bits)) {

          for (auto postfix : postfixes) {
            if (bitops::popcnt(postfix & mask) & 1) {
              vec_out(idx) += val_diff * vec_in(idx);
            } else {
              vec_out(idx) += val_same * vec_in(idx);
            }
            ++idx;
          }
        }

        // s2 is prefix s1 is postfix
        else {
          bit_t s1mask = (bit_t)1 << s1;
          bit_t s2mask = (bit_t)1 << (s2 - n_postfix_bits);

          // s2 is up
          if (prefix & s2mask) {
            for (auto postfix : postfixes) {
              if (postfix & s1mask) {
                vec_out(idx) += val_same * vec_in(idx);
              } else {
                vec_out(idx) += val_diff * vec_in(idx);
              }
              ++idx;
            }
          }

          // s2 is dn
          else {
            for (auto postfix : postfixes) {
              if (postfix & s1mask) {
		vec_out(idx) += val_diff * vec_in(idx);
              } else {
		vec_out(idx) += val_same * vec_in(idx);
              }
              ++idx;
            }
          }
        }
      }
    }
  }
}
} // namespace hydra::spinhalfterms
