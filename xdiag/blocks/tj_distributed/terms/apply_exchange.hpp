#pragma once

#include <xdiag/bits/bitops.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/parallel/mpi/buffer.hpp>
#include <xdiag/parallel/mpi/communicator.hpp>

#include <xdiag/blocks/tj_distributed/terms/generic_term_mixed.hpp>

namespace xdiag::tj_distributed {

template <typename bit_t, typename coeff_t, class Basis>
void apply_exchange(Op const &op, Basis &&basis, const coeff_t *vec_in,
                    coeff_t *vec_out) {
  using namespace bits;

  assert(op.coupling_defined());
  assert(op.type_defined());
  assert(op.size() == 2);
  assert(op.sites_disjoint());
  std::string type = op.type();
  assert(type == "EXCHANGE");

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  Coupling cpl = op.coupling();
  assert(cpl.isexplicit() && !cpl.ismatrix());
  coeff_t J = cpl.as<coeff_t>();
  coeff_t Jhalf = J / 2.;
  coeff_t Jhalf_conj = conj(Jhalf);

  int64_t s1 = op[0];
  int64_t s2 = op[1];

  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

  int64_t l = std::min(s1, s2);
  int64_t u = std::max(s1, s2);
  bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

  int64_t n_sites = basis.n_sites();
  int64_t n_up = basis.n_up();
  int64_t n_dn = basis.n_dn();
  int64_t n_dn_configurations = combinatorics::binomial(n_sites - n_up, n_dn);
  // Find out how many states is sent to each process
  std::vector<int64_t> n_states_i_send(mpi_size, 0);

  // Flip states and check out how much needs to be communicated
  int64_t idx_up = 0;
  for (bit_t up : basis.my_ups()) {
    if (popcnt(up & flipmask) == 1) {
      bit_t flipped_up = up ^ flipmask;
      int target = basis.rank(flipped_up);

      for (bit_t dn : basis.my_dns_for_ups(idx_up)) {
        if (popcnt(dn & flipmask) == 1)
          ++n_states_i_send[target];
      }
    }
    ++idx_up;
  }

  // Exchange information on who sends how much to whom
  mpi::Communicator comm(n_states_i_send);
  mpi::buffer.reserve<coeff_t>(comm.send_buffer_size(),
                               comm.recv_buffer_size());

  // Flip states and check how much needs to be communicated
  idx_up = 0;
  int64_t idx = 0;
  for (bit_t up : basis.my_ups()) {
    if (popcnt(up & flipmask) == 1) {
      bit_t flipped_up = up ^ flipmask;
      int target = basis.rank(flipped_up);

      for (bit_t dn : basis.my_dns_for_ups(idx_up)) {
        if (popcnt(dn & flipmask) == 1) {
          comm.add_to_send_buffer(target, vec_in[idx],
                                  mpi::buffer.send<coeff_t>());
        }
        ++idx;
      }
    } else {
      idx += n_dn_configurations;
    }
    ++idx_up;
  }

  // Alltoall called
  comm.all_to_all(mpi::buffer.send<coeff_t>(), mpi::buffer.recv<coeff_t>());

  // Get the original upspin configuration and its source proc
  std::vector<std::vector<bit_t>> ups_i_get_from_proc(mpi_size);
  for (bit_t up : basis.my_ups()) {
    bit_t flipped_up = up ^ flipmask;
    int source = basis.rank(flipped_up);
    ups_i_get_from_proc[source].push_back(up);
  }

  int64_t recv_idx = 0;

  // Loop over origin of arrived
  for (int m = 0; m < mpi_size; ++m) {

    // Sort according to order of flipped upspins
    std::sort(ups_i_get_from_proc[m].begin(), ups_i_get_from_proc[m].end(),
              [&flipmask](bit_t const &a, bit_t const &b) {
                bit_t flipped_a = a ^ flipmask;
                bit_t flipped_b = b ^ flipmask;
                return flipped_a < flipped_b;
              });

    for (bit_t up : ups_i_get_from_proc[m]) {
      if (popcnt(up & flipmask) == 1) {
        bool fermi_up = bits::popcnt(up & fermimask) & 1;
        bool up_s1_set = bits::gbit(up, s2);

        int64_t up_offset = basis.my_ups_offset(up);

        for (int64_t target_idx = up_offset;
             target_idx < up_offset + n_dn_configurations; ++target_idx) {
          bit_t dn = basis.my_dns_for_ups_storage(target_idx);

          if (bits::popcnt(dn & flipmask) == 1) {
            bool fermi_dn = bits::popcnt(dn & fermimask) & 1;

            if constexpr (isreal<coeff_t>()) {
              vec_out[target_idx] += ((fermi_up ^ fermi_dn) ? Jhalf : -Jhalf) *
                                     mpi::buffer.recv<coeff_t>()[recv_idx];
            } else {
              if (up_s1_set) {
                vec_out[target_idx] +=
                    ((fermi_up ^ fermi_dn) ? Jhalf : -Jhalf) *
                    mpi::buffer.recv<coeff_t>()[recv_idx];
              } else {
                vec_out[target_idx] +=
                    ((fermi_up ^ fermi_dn) ? Jhalf_conj : -Jhalf_conj) *
                    mpi::buffer.recv<coeff_t>()[recv_idx];
              }
            }
            ++recv_idx;
          }
        }
      }
    }

  } // loop over processes
}

} // namespace xdiag::tj_distributed
