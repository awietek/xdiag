#pragma once

#include <xdiag/bitops/bitops.hpp>
#include <xdiag/blocks/spinhalf_mpi/spinhalf_mpi.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/common.hpp>
#include <xdiag/mpi/logger_mpi.hpp>
#include <xdiag/operators/bondlist.hpp>
#include <xdiag/operators/couplings.hpp>

namespace xdiag::terms {

template <class bit_t, class coeff_t>
void spinhalf_mpi_ising(BondList const &bonds, Couplings const &couplings,
                        indexing::SpinhalfMPIIndexingSz<bit_t> const &indexing,
                        lila::Vector<coeff_t> const &vec_in,
                        lila::Vector<coeff_t> &vec_out) {

  auto clean_bonds =
      utils::clean_bondlist(bonds, couplings, {"HEISENBERG", "HB", "ISING"}, 2);

  for (auto bond : clean_bonds) {

    // Set values for same/diff
    std::string cpl = bond.coupling();
    coeff_t J = utils::get_coupling<coeff_t>(couplings, cpl);
    coeff_t val_same = J / 4.;
    coeff_t val_diff = -J / 4.;

    int s1 = std::min(bond.site(0), bond.site(1));
    int s2 = std::max(bond.site(0), bond.site(1));
    if (s1 == s2)
      LogMPI.err("Error computing SpinhalfMPI Ising: "
                 "operator acting on twice the same site");
    bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    int n_postfix_bits = indexing.n_postfix_bits();

    int64_t idx = 0;
    for (auto prefix : indexing.prefixes()) {

      bit_t prefix_shifted = (prefix << n_postfix_bits);
      auto postfixes = indexing.postfixes(prefix);

      // Both sites are on prefixes
      if ((s1 >= n_postfix_bits) && (s2 >= n_postfix_bits)) {
        coeff_t val =
            (bitops::popcnt(prefix_shifted & mask) & 1) ? val_diff : val_same;
        int64_t end = idx + postfixes.size();
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
} // namespace xdiag::terms
