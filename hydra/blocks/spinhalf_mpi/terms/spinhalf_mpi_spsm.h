#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/spinhalf_mpi/spinhalf_mpi.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/mpi/logger_mpi.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms {

template <class bit_t, class coeff_t>
void spinhalf_mpi_spsm(BondList const &bonds, Couplings const &couplings,
                       SpinhalfMPI<bit_t> const &block,
                       lila::Vector<coeff_t> const &vec_in,
                       lila::Vector<coeff_t> &vec_out, std::string spsm) {
  assert((spsm == "S+") || (spsm == "S-"));
  int n_up = block.n_up();
  int n_prefix_bits = block.n_prefix_bits_;
  int n_postfix_bits = block.n_postfix_bits_;

  auto clean_bonds = utils::clean_bondlist(bonds, couplings, {spsm}, 1);
  auto [prefix_bonds, postfix_bonds, mixed_bonds] =
    get_prefix_postfix_mixed_bonds(clean_bonds, n_postfix_bits);
  assert(mixed_bonds.size() == 0);

  // Apply S+S- if site belongs to postfixes
  for (auto bond : postfix_bonds) {

    int s = bond[0];
    assert(s < n_postfix_bits);
    bit_t mask = ((bit_t)1 << s);
    std::string cpl = bond.coupling();
    coeff_t H = utils::get_coupling<coeff_t>(couplings, cpl);

    idx_t idx = 0;
    idx_t idx_prefix = 0;
    for (auto prefix : block.prefixes_) {
      int n_up_prefix = bitops::popcnt(prefix);
      int n_up_postfix = block.n_up() - n_up_prefix;

      auto const &postfixes = block.postfix_states_[n_up_postfix];
      auto const &lintable = block.postfix_lintables_[n_up_postfix];

      if (spsm == "S+") {
        for (auto postfix_in : postfixes) {
          if (!(postfix_in & mask)) {
            bit_t postfix_out = postfix_in | mask;
            idx_t idx_out = idx_prefix + lintable.index(postfix_out);
            vec_out(idx_out) += H * vec_in(idx);
          }
          ++idx;
        }
      } else if (spsm == "S-") {
        for (auto postfix_in : postfixes) {
          if (postfix_in & mask) {
            bit_t postfix_out = postfix_in ^ mask;
            idx_t idx_out = idx_prefix + lintable.index(postfix_out);
            vec_out(idx_out) += H * vec_in(idx);
          }
          ++idx;
        }
      }
      idx_prefix += combinatorics::binomial(n_postfix_bits, n_up_postfix);
      assert(idx_prefix == idx);
    }
  }

  // Apply S+S- if site belongs to prefixes
  auto tpre = rightnow_mpi();
  if (prefix_bonds.size() > 0) {
    std::vector<coeff_t> send_buffer(vec_in.size(), 0);
    std::vector<coeff_t> recv_buffer(vec_in.size(), 0);

    // Transpose prefix/postfix order
    auto ttrans1 = rightnow_mpi();
    spinhalf_mpi_transpose(block, vec_in.vector(), send_buffer, recv_buffer,
                           false);
    timing_mpi(ttrans1, rightnow_mpi(), " (spsm) transpose 1", 2);

    // Loop over prefix_bonds
    for (auto bond : prefix_bonds) {
      assert(n_postfix_bits <= bond[0]);
      int s = bond[0] - n_postfix_bits;
      bit_t mask = ((bit_t)1 << s);
      std::string cpl = bond.coupling();
      coeff_t H = utils::get_coupling<coeff_t>(couplings, cpl);

      // loop through all postfixes
      idx_t idx = 0;
      idx_t idx_postfix = 0;
      for (auto postfix : block.postfixes_) {
        int n_up_postfix = bitops::popcnt(postfix);
        int n_up_prefix = n_up - n_up_postfix;

        auto const &lintable = block.prefix_lintables_[n_up_prefix];
        auto const &prefixes = block.prefix_states_[n_up_prefix];

        if (spsm == "S+") {
          for (auto prefix_in : prefixes) {
            if (!(prefix_in & mask)) {
              bit_t prefix_out = prefix_in | mask;
              idx_t idx_out = idx_postfix + lintable.index(prefix_out);
              recv_buffer[idx_out] += H * send_buffer[idx];
            }
            ++idx;
          }
        } else if (spsm == "S-") {
          for (auto prefix_in : prefixes) {
            if (prefix_in & mask) {
              bit_t prefix_out = prefix_in ^ mask;
              idx_t idx_out = idx_postfix + lintable.index(prefix_out);
              recv_buffer[idx_out] += H * send_buffer[idx];
            }
            ++idx;
          }
        }

        idx_postfix += combinatorics::binomial(n_prefix_bits, n_up_prefix);
        assert(idx_postfix == idx);
      }
    }

    // Transpose back
    auto ttrans2 = rightnow_mpi();
    spinhalf_mpi_transpose(block, recv_buffer, send_buffer, recv_buffer, true);
    timing_mpi(ttrans2, rightnow_mpi(), " (spsm) transpose 2", 2);

    for (idx_t idx = 0; idx < vec_out.size(); ++idx) {
      vec_out(idx) += send_buffer[idx];
    }
  }

  timing_mpi(tpre, rightnow_mpi(), " (spsm) prefix", 2);
}

} // namespace hydra::terms
