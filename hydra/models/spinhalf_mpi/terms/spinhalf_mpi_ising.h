#pragma once

#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/models/spinhalf_mpi/spinhalf_mpi.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>

namespace hydra::spinhalfterms {

template <class bit_t, class coeff_t>
void do_ising_mpi(BondList const &bonds, Couplings const &couplings,
                  SpinhalfMPI<bit_t> const &block,
		  lila::Vector<coeff_t> const& vec_in,
		  lila::Vector<coeff_t> & vec_out) {

  auto ising = bonds.bonds_of_type("HEISENBERG") +
               bonds.bonds_of_type("ISING") + bonds.bonds_of_type("HB");

  for (auto bond : ising) {

    if (bond.size() != 2)
      HydraLog.err("Error computing SpinhalfMPI Ising: "
                   "bond must have exactly two sites defined");

    std::string coupling = bond.coupling();
    if (couplings.defined(coupling) &&
        !lila::close(couplings[coupling], (complex)0.)) {

      double J = lila::real(couplings[coupling]);

      // Set values for same/diff (tJ model definition)
      std::string type = bond.type();
      double val_same = J / 4.;
      double val_diff = -J / 4.;

      int s1 = bond.site(0);
      int s2 = bond.site(1);
      if (s1 == s2)
        HydraLog.err("Error computing SpinhalfMPI Ising: "
                     "operator acting on twice the same site");
      bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

      int n_postfix_bits = block.n_postfix_bits_;

      idx_t idx = 0;
      for (auto prefix : block.prefixes_) {
        int n_up_prefix = utils::popcnt(prefix);
        int n_up_postfix = block.n_up() - n_up_prefix;
        for (auto postfix : Combinations<bit_t>(n_postfix_bits, n_up_postfix)) {
          bit_t state = (prefix << n_postfix_bits) | postfix;
          if (utils::popcnt(state & mask) & 1)
            vec_out(idx) += val_diff * vec_in(idx);
          else
            vec_out(idx) += val_same * vec_in(idx);
        }
        ++idx;
      }
    }
  }
}

} // namespace hydra::spinhalfterms
