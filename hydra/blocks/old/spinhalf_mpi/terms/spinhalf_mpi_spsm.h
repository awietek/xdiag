#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/spinhalf_mpi/spinhalf_mpi.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/mpi/buffer.h>
#include <hydra/mpi/logger_mpi.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms {

template <class bit_t, class coeff_t>
void spinhalf_mpi_spsm(
    BondList const &bonds, Couplings const &couplings,
    indexing::SpinhalfMPIIndexingSz<bit_t> const &indexing_in,
    lila::Vector<coeff_t> const &vec_in,
    indexing::SpinhalfMPIIndexingSz<bit_t> const &indexing_out,
    lila::Vector<coeff_t> &vec_out, std::string spsm) {
  assert(indexing_in.size() == vec_in.size());
  assert(indexing_out.size() == vec_out.size());
  assert(indexing_out.n_sites() == indexing_in.n_sites());

  int n_postfix_bits = indexing_in.n_postfix_bits();

  auto clean_bonds = utils::clean_bondlist(bonds, couplings, {spsm}, 1);
  if (clean_bonds.size() > 0) {
    if (spsm == "S+") {
      if (indexing_in.n_up() + 1 != indexing_out.n_up()) {
        LogMPI.err("Error: Invalid nup for S+! in: {}, out: {}",
                   indexing_in.n_up(), indexing_out.n_up());
      }
    } else if (spsm == "S-") {
      if (indexing_in.n_up() - 1 != indexing_out.n_up()) {
        LogMPI.err("Error: Invalid nup for S-! in: {}, out: {}",
                   indexing_in.n_up(), indexing_out.n_up());
      }
    }
    assert((spsm == "S+") || (spsm == "S-"));
  }

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
    for (auto prefix : indexing_in.prefixes()) {

      auto const &postfixes = indexing_in.postfixes(prefix);
      auto const &lintable = indexing_out.postfix_indexing(prefix);

      // lila::Log("nup: {} spsm: {}", n_up, spsm);
      // lila::Log("postfs: {}", postfixes.size());
      // lila::Log("prefix: {}", prefix);
      idx_t idx_prefix = indexing_out.prefix_begin(prefix);

      if (spsm == "S+") {
        for (auto postfix_in : postfixes) {
          if (!(postfix_in & mask)) {
            bit_t postfix_out = postfix_in | mask;
            idx_t idx_out = idx_prefix + lintable.index(postfix_out);
            // lila::Log("S+ s: {} in: {} {};{} -> out: {} {};{}", s, idx,
            //           bitops::bits_to_string(prefix, n_prefix_bits),
            //           bitops::bits_to_string(postfix_in, n_postfix_bits),
            //           idx_out, bitops::bits_to_string(prefix, n_prefix_bits),
            //           bitops::bits_to_string(postfix_out, n_postfix_bits));

            vec_out(idx_out) += H * vec_in(idx);
          }
          ++idx;
        }
        // idx_prefix += combinatorics::binomial(n_postfix_bits, n_up_postfix +
        // 1);
      } else if (spsm == "S-") {
        for (auto postfix_in : postfixes) {
          if (postfix_in & mask) {
            bit_t postfix_out = postfix_in ^ mask;
            idx_t idx_out = idx_prefix + lintable.index(postfix_out);
            // lila::Log("S- s: {} in: {} {};{} -> out: {} {};{}", s, idx,
            //           bitops::bits_to_string(prefix, n_prefix_bits),
            //           bitops::bits_to_string(postfix_in, n_postfix_bits),
            //           idx_out, bitops::bits_to_string(prefix, n_prefix_bits),
            //           bitops::bits_to_string(postfix_out, n_postfix_bits));

            // lila::Log("idx_out: {} = {} + {}", idx_out, idx_prefix,
            //           lintable.index(postfix_out));
            vec_out(idx_out) += H * vec_in(idx);
          }
          ++idx;
        }
        // lila::Log("comb: {}, {} {}",
        //           combinatorics::binomial(n_postfix_bits, n_up_postfix - 1),
        //           n_postfix_bits, n_up_postfix - 1);
        // idx_prefix += combinatorics::binomial(n_postfix_bits, n_up_postfix -
        // 1);
      }
    }
  }

  // Apply S+S- if site belongs to prefixes
  auto tpre = rightnow_mpi();
  if (prefix_bonds.size() > 0) {

    // Get communication buffers
    idx_t buffer_size =
        std::max(indexing_in.size_max(), indexing_out.size_max());
    mpi::buffer.reserve<coeff_t>(buffer_size);
    auto send_buffer = mpi::buffer.send<coeff_t>();
    auto recv_buffer = mpi::buffer.recv<coeff_t>();

    // Transpose prefix/postfix order
    auto ttrans1 = rightnow_mpi();
    spinhalf_mpi_transpose(indexing_in, vec_in.vector(), false);
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
      for (auto postfix : indexing_in.postfixes()) {

        auto const &prefixes = indexing_in.prefixes(postfix);
        auto const &lintable = indexing_out.prefix_indexing(postfix);
        idx_t idx_postfix = indexing_out.postfix_begin(postfix);

        if (spsm == "S+") {
          for (auto prefix_in : prefixes) {
            if (!(prefix_in & mask)) {
              bit_t prefix_out = prefix_in | mask;
              idx_t idx_out = idx_postfix + lintable.index(prefix_out);

              // int n_sites = indexing_in.n_sites();
              // LogMPI.out("s: {}, bond[0]: {}, n_postfix_bits: {}", s,
              // bond[0],
              //            n_postfix_bits);
              // LogMPI.out("vec_in.size: {}, vec_out.size(): {}, max: {}",
              //            vec_in.size(), vec_out.size(),
              //            std::max(vec_in.size(), vec_out.size()));

              // LogMPI.out("prefix_in: {}, prefix_out: {}, mask: {}",
              //            bitops::bits_to_string(prefix_in, n_prefix_bits),
              //            bitops::bits_to_string(prefix_out, n_prefix_bits),
              //            bitops::bits_to_string(mask, n_prefix_bits));

              // LogMPI.out("{} = {} (idx_postfix) + {} (index(prefix out)",
              //            idx_out, idx_postfix, lintable.index(prefix_out));

              // if (idx_out >= recv_buffer.size()) {
              //   int n_prefix_bits = indexing_in.n_prefix_bits();

              //   lila::Log.out(
              //       "recv_buffer[idx_out] {} {} send_buffer[idx] {} {}",
              //       idx_out, recv_buffer.size(), idx, send_buffer.size());
              //   lila::Log.out("s: {}, bond[0]: {}, n_postfix_bits: {}", s,
              //                 bond[0], n_postfix_bits);
              //   lila::Log.out("vec_in.size: {}, vec_out.size(): {}, max: {}",
              //                 vec_in.size(), vec_out.size(),
              //                 std::max(vec_in.size(), vec_out.size()));
              //   lila::Log.out("postfix: {}",
              //                 bitops::bits_to_string(postfix,
              //                 n_postfix_bits));

              //   lila::Log.out("prefix_in: {}, prefix_out: {}, mask: {}",
              //                 bitops::bits_to_string(prefix_in,
              //                 n_prefix_bits),
              //                 bitops::bits_to_string(prefix_out,
              //                 n_prefix_bits), bitops::bits_to_string(mask,
              //                 n_prefix_bits));

              //   lila::Log("my postfixes");
              //   for (auto pf : indexing_out.postfixes()) {
              //     lila::Log("pf: {}  begin: {} size: {}",
              //               bitops::bits_to_string(pf, n_postfix_bits),
              //               indexing_out.postfix_begin(pf),
              //               indexing_out.size());
              //   }

              //   lila::Log.out(
              //       "#{} {} = {} (idx_postfix) + {} (index(prefix out)",
              //       indexing_out.mpi_rank(), idx_out, idx_postfix,
              //       lintable.index(prefix_out));
              //   assert(idx_out < recv_buffer.size());
              // }
              recv_buffer[idx_out] += H * send_buffer[idx];
            }
            ++idx;
          }
          // idx_postfix +=
          //     combinatorics::binomial(n_prefix_bits, n_up_prefix + 1);
        } else if (spsm == "S-") {
          for (auto prefix_in : prefixes) {
            if (prefix_in & mask) {
              bit_t prefix_out = prefix_in ^ mask;
              idx_t idx_out = idx_postfix + lintable.index(prefix_out);
              recv_buffer[idx_out] += H * send_buffer[idx];
            }
            ++idx;
          }
          // idx_postfix +=
          //     combinatorics::binomial(n_prefix_bits, n_up_prefix - 1);
        }
      }
    }

    // Transpose back
    auto ttrans2 = rightnow_mpi();
    spinhalf_mpi_transpose(indexing_out, recv_buffer, true);
    timing_mpi(ttrans2, rightnow_mpi(), " (spsm) transpose 2", 2);

    for (idx_t idx = 0; idx < vec_out.size(); ++idx) {
      vec_out(idx) += send_buffer[idx];
    }
  }

  timing_mpi(tpre, rightnow_mpi(), " (spsm) prefix", 2);
}

} // namespace hydra::terms
