#include "tjmodelmpi.h"

#include "hubbardmodeldetail.h"
#include <hydra/utils/bitops.h>
#include <hydra/combinatorics/up_down_hole.h>
#include <hydra/combinatorics/binomial.h>

namespace hydra {

using namespace combinatorics;

template <class coeff_t, class bit_t, class idx_t>
TJModelMPI<coeff_t, bit_t, idx_t>::TJModelMPI(BondList bondlist,
                                              Couplings couplings, qn_tj qn)
    : n_sites_(bondlist.n_sites()), qn_(qn) {
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

  // Currently unused operatuors
  std::vector<std::pair<int, int>> currents_;
  std::vector<coeff_t> current_amplitudes_;
  std::vector<std::pair<int, int>> interactions_;
  std::vector<double> interaction_strengths_;
  double U_;

  // Use Hubbard routine so set interation terms
  detail::set_hubbard_terms<coeff_t>(
      bondlist, couplings, hoppings_, hopping_amplitudes_, currents_,
      current_amplitudes_, interactions_, interaction_strengths_, onsites_,
      onsite_potentials_, szszs_, szsz_amplitudes_, exchanges_,
      exchange_amplitudes_, U_);

  initialize();
}

template <class coeff_t, class bit_t, class idx_t>
void TJModelMPI<coeff_t, bit_t, idx_t>::apply_hamiltonian(
    const lila::VectorMPI<coeff_t> &in_vec, lila::VectorMPI<coeff_t> &out_vec,
    bool verbose) {
  using utils::gbit;
  using utils::gbits;
  using utils::popcnt;
  using namespace hydra::combinatorics;

  if ((mpi_rank_ == 0) && verbose) {
    printf("max_local_dim: %lu\n", max_local_dim_);
    printf("min_local_dim: %lu\n", min_local_dim_);
  }

  Zeros(out_vec.vector_local());

  // Apply szsz terms
  double t1 = MPI_Wtime();
  int szsz_idx = 0;
  for (auto pair : szszs_) {
    int s1 = std::min(pair.first, pair.second);
    int s2 = std::max(pair.first, pair.second);
    double jz = szsz_amplitudes_[szsz_idx] * 0.25;
    bit_t s1_mask = (bit_t)1 << s1;
    bit_t s2_mask = (bit_t)1 << s2;

    if (std::abs(jz) > 1e-14) {
      idx_t n_hole_configurations = hs_holes_in_ups_.size();
      idx_t upspin_idx = 0;
      for (bit_t upspins : my_upspins_) {
        idx_t upspin_offset = upspin_idx * n_hole_configurations;
        bit_t upspin_s1 = upspins & s1_mask;
        bit_t upspin_s2 = upspins & s2_mask;

        // Both upspins are set SzSz -> +0.25
        if (upspin_s1 && upspin_s2) {
          for (idx_t idx = upspin_offset;
               idx < upspin_offset + n_hole_configurations; ++idx) {
            out_vec(idx) += jz * in_vec(idx);
          }
        }
        // upspin at s1 is set -> -0.25 if downspin s2 is set
        else if (upspin_s1) {
          for (idx_t idx = upspin_offset;
               idx < upspin_offset + n_hole_configurations; ++idx) {
            bit_t downspins = downspins_table_[idx];
            out_vec(idx) -= jz * ((downspins & s2_mask) != 0) * in_vec(idx);
          }
        }
        // upspin at s2 is set -> -0.25 if downspin s1 is set
        else if (upspin_s2) {
          for (idx_t idx = upspin_offset;
               idx < upspin_offset + n_hole_configurations; ++idx) {
            bit_t downspins = downspins_table_[idx];
            out_vec(idx) -= jz * ((downspins & s1_mask) != 0) * in_vec(idx);
          }
        }
        // no upspin set -> +0.25 if both downspins are set
        else {
          for (idx_t idx = upspin_offset;
               idx < upspin_offset + n_hole_configurations; ++idx) {
            bit_t downspins = downspins_table_[idx];
            out_vec(idx) += jz *
                            ((downspins & s1_mask) && (downspins & s2_mask)) *
                            in_vec(idx);
          }
        }

        ++upspin_idx;

      } // (bit_t upspins : my_upspins_)
    }
    ++szsz_idx;
  }
  double t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  diag: %3.4f\n", t2 - t1);

  // Spin exchange terms
  int exchange_idx = 0;
  t1 = MPI_Wtime();
  for (auto pair : exchanges_) {
    const int s1 = std::min(pair.first, pair.second);
    const int s2 = std::max(pair.first, pair.second);
    const coeff_t jx = exchange_amplitudes_[exchange_idx] * 0.5;
    const bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

    idx_t n_hole_configurations = hs_holes_in_ups_.size();

    // Find out how many states is sent to each process
    std::vector<int> n_states_i_send(mpi_size_, 0);
    if (std::abs(jx) > 1e-14) {
      // Flip states and check out how much needs to be communicated
      idx_t upspin_idx = 0;
      for (bit_t upspins : my_upspins_) {
        if (popcnt(upspins & flipmask) == 1) {
          idx_t upspin_offset = upspin_idx * n_hole_configurations;

          // idx_t upspin_offset2 = my_upspins_offset_[upspins];
          // assert(upspin_offset == upspin_offset2);

          bit_t flipped_upspins = upspins ^ flipmask;
          int target = mpi_rank_of_spins(flipped_upspins);

          for (idx_t idx = upspin_offset;
               idx < upspin_offset + n_hole_configurations; ++idx) {
            bit_t downspins = downspins_table_[idx];
            if (popcnt(downspins & flipmask) == 1)
              ++n_states_i_send[target];
          }
        }
        ++upspin_idx;
      }

      // Exchange information on who sends how much to whom
      std::vector<int> n_states_i_recv(mpi_size_, 0);
      lila::MPI_Alltoall<int>(n_states_i_send.data(), 1, n_states_i_recv.data(),
                              1, MPI_COMM_WORLD);

      // Sum up states sent/recvd and eventually resize buffers
      idx_t sum_n_states_i_send = std::accumulate(
          n_states_i_send.begin(), n_states_i_send.end(), (idx_t)0);
      idx_t sum_n_states_i_recv = std::accumulate(
          n_states_i_recv.begin(), n_states_i_recv.end(), (idx_t)0);

      if (sum_n_states_i_send > send_buffer_.size())
        send_buffer_.resize(sum_n_states_i_send);
      if (sum_n_states_i_recv > recv_buffer_.size())
        recv_buffer_.resize(sum_n_states_i_recv);

      // Compute offsets of send/recieve states
      std::vector<int> n_states_i_send_offsets(mpi_size_, 0);
      std::vector<int> n_states_i_recv_offsets(mpi_size_, 0);
      for (int m = 0; m < mpi_size_; ++m)
        for (int n = 0; n < m; ++n) {
          n_states_i_send_offsets[m] += n_states_i_send[n];
          n_states_i_recv_offsets[m] += n_states_i_recv[n];
        }

      // Flip states and check how much needs to be communicated
      std::vector<int> n_states_prepared(mpi_size_, 0);
      upspin_idx = 0;
      for (bit_t upspins : my_upspins_) {
        if (popcnt(upspins & flipmask) == 1) {
          idx_t upspin_offset = upspin_idx * n_hole_configurations;

          // idx_t upspin_offset2 = my_upspins_offset_[upspins];
          // assert(upspin_offset == upspin_offset2);

          bit_t flipped_upspins = upspins ^ flipmask;
          int target = mpi_rank_of_spins(flipped_upspins);

          for (idx_t idx = upspin_offset;
               idx < upspin_offset + n_hole_configurations; ++idx) {
            bit_t downspins = downspins_table_[idx];
            if (popcnt(downspins & flipmask) == 1) {
              int send_idx =
                  n_states_i_send_offsets[target] + n_states_prepared[target];
              send_buffer_[send_idx] = in_vec.vector_local()(idx);
              ++n_states_prepared[target];
            }
          }
        }
        ++upspin_idx;
      }

      // Check, whether correct number of states has been prepared
      for (int m = 0; m < mpi_size_; ++m)
        assert(n_states_prepared[m] == n_states_i_send[m]);

      // Alltoall call
      lila::MPI_Alltoallv<coeff_t>(send_buffer_.data(), n_states_i_send.data(),
                                   n_states_i_send_offsets.data(),
                                   recv_buffer_.data(), n_states_i_recv.data(),
                                   n_states_i_recv_offsets.data(),
                                   MPI_COMM_WORLD);

      // Get the original upspin configuration and its source proc
      std::vector<std::vector<bit_t>> upspins_i_get_from_proc(mpi_size_);
      for (const bit_t &upspins : my_upspins_) {
        bit_t flipped_upspins = upspins ^ flipmask;
        int source = mpi_rank_of_spins(flipped_upspins);
        upspins_i_get_from_proc[source].push_back(upspins);
      }

      idx_t recv_idx = 0;
      for (int m = 0; m < mpi_size_; ++m) {

        // Sort according to order of flipped upspins
        std::sort(upspins_i_get_from_proc[m].begin(),
                  upspins_i_get_from_proc[m].end(),
                  [&flipmask](bit_t const &a, bit_t const &b) {
                    bit_t flipped_a = a ^ flipmask;
                    bit_t flipped_b = b ^ flipmask;
                    return flipped_a < flipped_b;
                  });

        for (const bit_t &upspins : upspins_i_get_from_proc[m]) {
          if (popcnt(upspins & flipmask) == 1) {
            idx_t upspin_offset = my_upspins_offset_[upspins];
            for (idx_t target_idx = upspin_offset;
                 target_idx < upspin_offset + n_hole_configurations;
                 ++target_idx) {
              auto downspins = downspins_table_[target_idx];
              if (popcnt(downspins & flipmask) == 1) {
                out_vec.vector_local()(target_idx) +=
                    jx * recv_buffer_[recv_idx];
                ++recv_idx;
              }
            }
          }
        }

      } // loop over processes
    }   // if (std::abs(jx) > 1e-14)
    ++exchange_idx;
  } // for (auto pair: exchanges_)
  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  exch: %3.4f\n", t2 - t1);

  // Apply hoppings on downspins
  t1 = MPI_Wtime();
  int hopping_idx = 0;
  for (auto pair : hoppings_) {
    int s1 = std::min(pair.first, pair.second);
    int s2 = std::max(pair.first, pair.second);
    coeff_t t = hopping_amplitudes_[hopping_idx];
    if (pair.first > pair.second) t = lila::conj(t);
    
    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    bit_t firstmask = (bit_t)1 << s1;

    idx_t n_hole_configurations = hs_holes_in_ups_.size();

    if (std::abs(t) > 1e-14) {
      // Loop over all configurations
      idx_t upspin_idx = 0;
      for (bit_t upspins : my_upspins_) {
        // t-J hard-core constraint
        if ((upspins & flipmask) == 0) {
          idx_t upspins_offset = upspin_idx * n_hole_configurations;

          for (idx_t idx = upspins_offset;
               idx < upspins_offset + n_hole_configurations; ++idx) {
            bit_t downspins = downspins_table_[idx];

            // Check if hopping is possible
            if (((downspins & flipmask) != 0) &&
                ((downspins & flipmask) != flipmask)) {
              double fermi = popcnt(gbits(upspins ^ downspins, s2 - s1, s1)) & 1
                                 ? -1.
                                 : 1.;
              bit_t new_downspins = downspins ^ flipmask;

              // Binary search the downspin configuration
              auto it =
                  std::lower_bound(downspins_table_.begin() + upspins_offset,
                                   downspins_table_.begin() + upspins_offset +
                                       n_hole_configurations,
                                   new_downspins);
              idx_t new_idx = std::distance(downspins_table_.begin(), it);

              if (downspins & firstmask) {
                out_vec.vector_local()(new_idx) +=
                    fermi * t * in_vec.vector_local()(idx);
              } else {
                out_vec.vector_local()(new_idx) -=
                    fermi * lila::conj(t) * in_vec.vector_local()(idx);
              }
            }
          }
        } // if ((upspins & flipmask) == 0)

        ++upspin_idx;

      } // for(const bit_t& upspins : my_upspins_)
    }   // if (std::abs(t) > 1e-14)
    ++hopping_idx;
  } // hopping on downspin
  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  down: %3.4f\n", t2 - t1);

  // Send configurations
  t1 = MPI_Wtime();
  idx_t n_hole_configurations = hs_holes_in_ups_.size();
  std::vector<idx_t> n_states_already_prepared(mpi_size_, 0);
  for (idx_t down_idx = 0; down_idx < n_hole_configurations; ++down_idx) {
    for (idx_t up_idx = 0; up_idx < my_upspins_.size(); ++up_idx) {
      // bit_t upspins = my_upspins_[up_idx];
      // idx_t upspin_offset = my_upspins_offset_[upspins];

      idx_t upspin_offset = up_idx * n_hole_configurations;
      idx_t idx = upspin_offset + down_idx;

      bit_t downspins = downspins_table_[idx];
      int destination_mpi_rank = mpi_rank_of_spins(downspins);

      idx_t send_idx =
          n_downspins_i_send_forward_offsets_[destination_mpi_rank] +
          n_states_already_prepared[destination_mpi_rank]++;
      send_buffer_[send_idx] = in_vec.vector_local()(idx);
    }
  }
  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  send forward (prepare): %3.4f\n", t2 - t1);

  t1 = MPI_Wtime();
  lila::MPI_Alltoallv<coeff_t>(
      send_buffer_.data(), n_downspins_i_send_forward_.data(),
      n_downspins_i_send_forward_offsets_.data(), recv_buffer_.data(),
      n_downspins_i_recv_forward_.data(),
      n_downspins_i_recv_forward_offsets_.data(), MPI_COMM_WORLD);
  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  send forward (Alltoallv): %3.4f\n", t2 - t1);

  // // DEBUG START
  // // Fat DEBUG Print
  // MPI_Barrier(MPI_COMM_WORLD);

  // for(int nt = 0; nt < mpi_size_; ++nt)
  //   {
  //     if(mpi_rank_ == nt)
  //       {
  // 	printf("[%d]:\n", mpi_rank_);

  // 	idx_t upspin_idx = 0;
  // 	for (bit_t upspins : my_upspins_)
  // 	  {
  // 	    idx_t upspin_offset = my_upspins_offset_[upspins];
  // 	    idx_t downspin_offset = 0;
  // 	    for (bit_t holes : hs_holes_in_ups_)
  // 	      {
  // 		bit_t downspins = up_hole_to_down(upspins, holes);
  // 		idx_t idx = upspin_offset + downspin_offset;
  // 		printf("%s; %s %f -> [%d]\n",
  // 		       PrintSpinhalf(n_sites_, upspins).c_str(),
  // 		       PrintSpinhalf(n_sites_, downspins).c_str(),
  // 		       in_vec.vector_local()(idx),
  // 		       mpi_rank_of_spins(downspins));
  // 		++downspin_offset;
  // 	      }
  // 	    ++upspin_idx;
  // 	  }

  // 	printf("\n");
  //       }
  //     MPI_Barrier(MPI_COMM_WORLD);
  //   }
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(mpi_rank_ == 0) printf("SEND BUFFER
  // --------------------------------\n"); MPI_Barrier(MPI_COMM_WORLD); for(int
  // nt = 0; nt < mpi_size_; ++nt)
  //   {
  //     if(mpi_rank_ == nt)
  //       {
  // 	      printf("[%d]:\n", mpi_rank_);

  // 	      for (int sendto = 0; sendto < mpi_size_; ++sendto)
  // 		{
  // 		  int from = n_downspins_i_send_forward_offsets_[sendto];
  // 		  int to = (sendto == mpi_size_-1) ?
  // sum_n_downspins_i_send_forward_ :
  // n_downspins_i_send_forward_offsets_[sendto+1]; 		  for (int idx =
  // from; idx < to; ++idx) 		    printf("%f -> [%d] %s; %s\n",
  // send_buffer_[idx], sendto, PrintSpinhalf(n_sites_,
  // upspins_i_send_forward_[idx]).c_str(),
  // PrintSpinhalf(n_sites_, downspins_i_send_forward_[idx]).c_str());
  // 		}
  // 	      printf("\n");
  //       }
  //     MPI_Barrier(MPI_COMM_WORLD);
  //   }
  // MPI_Barrier(MPI_COMM_WORLD);

  // if(mpi_rank_ == 0) printf("RECV BUFFER
  // --------------------------------\n"); MPI_Barrier(MPI_COMM_WORLD); for(int
  // nt = 0; nt < mpi_size_; ++nt)
  //   {
  //     if(mpi_rank_ == nt)
  //       {
  // 	printf("[%d]:\n", mpi_rank_);
  // 	for (int recvfrom = 0; recvfrom < mpi_size_; ++recvfrom)
  // 	  {
  // 	    int from = n_downspins_i_recv_forward_offsets_[recvfrom];
  // 	    int to = (recvfrom == mpi_size_-1) ? sum_n_downspins_i_recv_forward_
  // : n_downspins_i_recv_forward_offsets_[recvfrom+1]; 	    for (int idx
  // = from; idx < to; ++idx) 	      printf("%f <- [%d]  %s; %s\n",
  // recv_buffer_[idx], recvfrom, 		     PrintSpinhalf(n_sites_,
  // upspins_i_recv_forward_[idx]).c_str(), PrintSpinhalf(n_sites_,
  // downspins_i_recv_forward_[idx]).c_str());

  // 	  }
  // 	printf("\n");
  //       }
  //     MPI_Barrier(MPI_COMM_WORLD);
  //   }
  // MPI_Barrier(MPI_COMM_WORLD);
  // // DEBUG END

  // Sort the received data into the send buffer
  t1 = MPI_Wtime();
  std::fill(send_buffer_.begin(), send_buffer_.end(), 0);
  n_hole_configurations = hs_holes_in_downs_.size();
  for (idx_t idx = 0; idx < sum_n_downspins_i_recv_forward_; ++idx) {
    bit_t upspins = upspins_i_recv_forward_[idx];
    bit_t downspins = downspins_i_recv_forward_[idx];

    // Binary search the upspin configuration
    idx_t downspins_offset = my_downspins_offset_[downspins];
    auto it = std::lower_bound(upspins_table_.begin() + downspins_offset,
                               upspins_table_.begin() + downspins_offset +
                                   n_hole_configurations,
                               upspins);
    idx_t sorted_idx = std::distance(upspins_table_.begin(), it);

    // bit_t holes = down_up_to_hole(downspins, upspins);
    // idx_t sorted_idx = my_downspins_offset_[downspins] +
    //   indexing_holes_in_downs_.index(holes);
    send_buffer_[sorted_idx] = recv_buffer_[idx];
  }
  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  sort forward: %3.4f\n", t2 - t1);

  // // DEBUG START
  // if(mpi_rank_ == 0) printf("SORTED SEND BUFFER
  // --------------------------------\n"); MPI_Barrier(MPI_COMM_WORLD); for(int
  // nt = 0; nt < mpi_size_; ++nt)
  //   {
  //     if(mpi_rank_ == nt)
  //       {
  // 	printf("[%d]:\n", mpi_rank_);
  // 	idx_t idx = 0;
  // 	for (bit_t downspins : my_downspins_)
  // 	  {
  // 	    for (bit_t holes : hs_holes_in_downs_)
  // 	      {
  // 		bit_t upspins = down_hole_to_up(downspins, holes);
  // 		printf("%s; %s  %f\n",
  // 		       PrintSpinhalf(n_sites_, upspins).c_str(),
  // 		       PrintSpinhalf(n_sites_, downspins).c_str(),
  // 		       send_buffer_[idx]);
  // 		++idx;
  // 	      }
  // 	  }
  // 	printf("\n");
  //       }
  //     MPI_Barrier(MPI_COMM_WORLD);
  //   }
  // MPI_Barrier(MPI_COMM_WORLD);
  // // DEBUG END

  // // DEBUG START
  // // JUST COPY INVEC FOR DEBUGGING
  // std::copy(send_buffer_.begin(), send_buffer_.end(), recv_buffer_.begin());
  // // DEBUG END

  // Apply hoppings on upspins
  t1 = MPI_Wtime();
  hopping_idx = 0;
  std::fill(recv_buffer_.begin(), recv_buffer_.end(),
            0); // clear the recv buffer
  for (auto pair : hoppings_) {
    int s1 = std::min(pair.first, pair.second);
    int s2 = std::max(pair.first, pair.second);
    coeff_t t = hopping_amplitudes_[hopping_idx];
    if (pair.first > pair.second) t = lila::conj(t);

    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    bit_t firstmask = (bit_t)1 << s1;

    if (std::abs(t) > 1e-14) {
      // Loop over all configurations
      idx_t downspin_idx = 0;
      for (bit_t downspins : my_downspins_) {
        // t-J hard-core constraint
        if ((downspins & flipmask) == 0) {
          idx_t downspin_offset = downspin_idx * n_hole_configurations;
          for (idx_t idx = downspin_offset;
               idx < downspin_offset + n_hole_configurations; ++idx) {
            bit_t upspins = upspins_table_[idx];

            // Check if hopping is possible
            if (((upspins & flipmask) != 0) &&
                ((upspins & flipmask) != flipmask)) {
              double fermi = popcnt(gbits(upspins ^ downspins, s2 - s1, s1)) & 1
                                 ? -1.
                                 : 1.;
              bit_t new_upspins = upspins ^ flipmask;

              // Binary search the upspin configuration
              auto it =
                  std::lower_bound(upspins_table_.begin() + downspin_offset,
                                   upspins_table_.begin() + downspin_offset +
                                       n_hole_configurations,
                                   new_upspins);
              idx_t new_idx = std::distance(upspins_table_.begin(), it);

              if (upspins & firstmask) {
                recv_buffer_[new_idx] += fermi * t * send_buffer_[idx];
              } else {
                recv_buffer_[new_idx] -=
                    fermi * lila::conj(t) * send_buffer_[idx];
              }
            }
          }
        } // if ((downspins & flipmask) == 0)

        ++downspin_idx;

      } // for (bit_t downspins : my_downspins_)
    }   // if (std::abs(t) > 1e-14)

    ++hopping_idx;
  } // hopping on upspins
  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  up:   %3.4f\n", t2 - t1);

  // // DEBUG START
  // if(mpi_rank_ == 0) printf("RESULTING RECV BUFFER
  // --------------------------------\n"); MPI_Barrier(MPI_COMM_WORLD); for(int
  // nt = 0; nt < mpi_size_; ++nt)
  //   {
  //     if(mpi_rank_ == nt)
  //       {
  // 	printf("[%d]:\n", mpi_rank_);
  // 	idx_t idx = 0;
  // 	for (bit_t downspins : my_downspins_)
  // 	  {
  // 	    for (bit_t holes : hs_holes_in_downs_)
  // 	      {
  // 		bit_t upspins = down_hole_to_up(downspins, holes);
  // 		printf("%s; %s  %f -> [%d]\n",
  // 		       PrintSpinhalf(n_sites_, upspins).c_str(),
  // 		       PrintSpinhalf(n_sites_, downspins).c_str(),
  // 		       recv_buffer_[idx], mpi_rank_of_spins(upspins));
  // 		++idx;
  // 	      }
  // 	  }
  // 	printf("\n");
  //       }
  //     MPI_Barrier(MPI_COMM_WORLD);
  //   }
  // MPI_Barrier(MPI_COMM_WORLD);
  // // DEBUG END

  // Send back the resulting configurations
  t1 = MPI_Wtime();
  std::fill(n_states_already_prepared.begin(), n_states_already_prepared.end(),
            0);
  for (idx_t up_idx = 0; up_idx < n_hole_configurations; ++up_idx) {
    for (idx_t down_idx = 0; down_idx < my_downspins_.size(); ++down_idx) {
      idx_t downspin_offset = down_idx * n_hole_configurations;

      // bit_t downspins = my_downspins_[down_idx];
      // idx_t downspin_offset2 = my_downspins_offset_[downspins];
      // assert(downspin_offset == downspin_offset2);

      idx_t idx = downspin_offset + up_idx;

      bit_t upspins = upspins_table_[idx];
      int destination_mpi_rank = mpi_rank_of_spins(upspins);

      idx_t send_idx = n_upspins_i_send_back_offsets_[destination_mpi_rank] +
                       n_states_already_prepared[destination_mpi_rank]++;
      send_buffer_[send_idx] = recv_buffer_[idx];
    }
  }

  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  send back (prepare):   %3.4f\n", t2 - t1);

  t1 = MPI_Wtime();
  lila::MPI_Alltoallv<coeff_t>(
      send_buffer_.data(), n_upspins_i_send_back_.data(),
      n_upspins_i_send_back_offsets_.data(), recv_buffer_.data(),
      n_upspins_i_recv_back_.data(), n_upspins_i_recv_back_offsets_.data(),
      MPI_COMM_WORLD);
  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  send back (Alltoall):   %3.4f\n", t2 - t1);

  // // DEBUG START
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(mpi_rank_ == 0) printf("BACK RECV BUFFER
  // --------------------------------\n"); MPI_Barrier(MPI_COMM_WORLD); for(int
  // nt = 0; nt < mpi_size_; ++nt)
  //   {
  //     if(mpi_rank_ == nt)
  //       {
  // 	printf("[%d]:\n", mpi_rank_);
  // 	for (int recvfrom = 0; recvfrom < mpi_size_; ++recvfrom)
  // 	  {
  // 	    int from = n_upspins_i_recv_back_offsets_[recvfrom];
  // 	    int to = (recvfrom == mpi_size_-1) ? sum_n_upspins_i_recv_back_ :
  // n_upspins_i_recv_back_offsets_[recvfrom+1]; 	    for (int idx = from;
  // idx < to;
  // ++idx) 	      printf("%f <- [%d]  %s; %s\n", recv_buffer_[idx],
  // recvfrom, 		     PrintSpinhalf(n_sites_,
  // upspins_i_recv_back_[idx]).c_str(), PrintSpinhalf(n_sites_,
  // downspins_i_recv_back_[idx]).c_str());

  // 	  }
  // 	printf("\n");
  //       }
  //     MPI_Barrier(MPI_COMM_WORLD);
  //   }
  // MPI_Barrier(MPI_COMM_WORLD);
  // // DEBUG END

  // Sort the received data into the send buffer
  t1 = MPI_Wtime();
  std::fill(send_buffer_.begin(), send_buffer_.end(), 0);
  n_hole_configurations = hs_holes_in_ups_.size();
  for (idx_t idx = 0; idx < sum_n_upspins_i_recv_back_; ++idx) {
    bit_t upspins = upspins_i_recv_back_[idx];
    bit_t downspins = downspins_i_recv_back_[idx];

    // Binary search the downspin configuration
    idx_t upspins_offset = my_upspins_offset_[upspins];
    auto it = std::lower_bound(downspins_table_.begin() + upspins_offset,
                               downspins_table_.begin() + upspins_offset +
                                   n_hole_configurations,
                               downspins);
    idx_t sorted_idx = std::distance(downspins_table_.begin(), it);

    // bit_t holes = up_down_to_hole(upspins, downspins);
    // idx_t sorted_idx = my_upspins_offset_[upspins] +
    //   indexing_holes_in_ups_.index(holes);

    send_buffer_[sorted_idx] = recv_buffer_[idx];
  }

  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  sort back:   %3.4f\n", t2 - t1);

  // // DEBUG START
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(mpi_rank_ == 0) printf("BACK SORTED SEND BUFFER
  // --------------------------------\n"); MPI_Barrier(MPI_COMM_WORLD); for(int
  // nt = 0; nt < mpi_size_; ++nt)
  //   {
  //     if(mpi_rank_ == nt)
  //       {
  // 	printf("[%d]:\n", mpi_rank_);
  // 	idx_t idx = 0;
  // 	for (bit_t upspins : my_upspins_)
  // 	  {
  // 	    for (bit_t holes : hs_holes_in_ups_)
  // 	      {
  // 		bit_t downspins = up_hole_to_down(upspins, holes);
  // 		printf("%s; %s  %f\n",
  // 		       PrintSpinhalf(n_sites_, upspins).c_str(),
  // 		       PrintSpinhalf(n_sites_, downspins).c_str(),
  // 		       send_buffer_[idx]);
  // 		++idx;
  // 	      }
  // 	  }
  // 	printf("\n");
  //       }
  //     MPI_Barrier(MPI_COMM_WORLD);
  //   }
  // MPI_Barrier(MPI_COMM_WORLD);
  // // DEBUG END

  // // DEBUG START
  // // DEBUG FORWARD_BACKWARD COMMUNICATION
  // for (int k=0; k<in_vec.vector_local().size(); ++k)
  //   assert(send_buffer_[k] == in_vec.vector_local()(k));
  // printf("[%d] SUCCESS !!!!!!!!!!!!!!!!!\n", mpi_rank_);
  // MPI_Finalize();
  // exit(1);
  // // DEBUG END

  // Fill in from send buffer
  t1 = MPI_Wtime();
  for (int k = 0; k < in_vec.vector_local().size(); ++k)
    out_vec.vector_local()(k) += send_buffer_[k];
  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  fill outvec:   %3.4f\n", t2 - t1);
}

template <class coeff_t, class bit_t, class idx_t>
void TJModelMPI<coeff_t, bit_t, idx_t>::apply_sz(
    const lila::VectorMPI<coeff_t> &in_vec, lila::VectorMPI<coeff_t> &out_vec,
    int site) {
  using utils::gbit;
  using utils::gbits;
  using utils::popcnt;
  assert(site < n_sites_);

  // Allocate out_vec
  try {
    out_vec.resize(local_dim_);
  } catch (...) {
    std::cerr << "[ " << mpi_rank_ << " ] Error: Could not allocate out_vector!"
              << std::endl
              << std::flush;
    MPI_Abort(MPI_COMM_WORLD, 4);
  }

  bit_t spin_mask = (bit_t)1 << site;

  idx_t n_hole_configurations = hs_holes_in_ups_.size();
  idx_t upspin_idx = 0;
  for (bit_t upspins : my_upspins_) {
    idx_t upspin_offset = upspin_idx * n_hole_configurations;
    bit_t upspin = upspins & spin_mask;

    // Upspin is set -> *0.5
    if (upspin) {
      for (idx_t idx = upspin_offset;
           idx < upspin_offset + n_hole_configurations; ++idx) {
        out_vec(idx) += 0.5 * in_vec(idx);
      }
      // Upspin isn't set -> iterate through downspins, *-0.5 if downspin is set
    } else {
      for (idx_t idx = upspin_offset;
           idx < upspin_offset + n_hole_configurations; ++idx) {
        bit_t downspins = downspins_table_[idx];
        out_vec(idx) -= 0.5 * ((downspins & spin_mask) != 0) * in_vec(idx);
      }
    }
    ++upspin_idx;
  }
}

template <class coeff_t, class bit_t, class idx_t>
void TJModelMPI<coeff_t, bit_t, idx_t>::initialize() {
  assert(valid(qn_, n_sites_));

  hs_upspins_ = BasisSpinHalf<bit_t>(n_sites_, qn_.n_up);
  hs_downspins_ = BasisSpinHalf<bit_t>(n_sites_, qn_.n_dn);
  hs_holes_in_ups_ = BasisSpinHalf<bit_t>(n_sites_ - qn_.n_up, qn_.n_dn);
  hs_holes_in_downs_ = BasisSpinHalf<bit_t>(n_sites_ - qn_.n_dn, qn_.n_up);

  indexing_holes_in_ups_ =
      LinTable<BasisSpinHalf<bit_t>, idx_t>(hs_holes_in_ups_);
  indexing_holes_in_downs_ =
      LinTable<BasisSpinHalf<bit_t>, idx_t>(hs_holes_in_downs_);

  // clear previous upspin, downspin members
  my_upspins_.clear();
  my_upspins_offset_.clear();
  my_downspins_.clear();
  my_downspins_offset_.clear();

  // Collect all upspin configurations belonging to this mpi_rank
  idx_t offset = 0;
  for (auto state : hs_upspins_) {
    if (mpi_rank_of_spins(state.spins) == mpi_rank_) {
      my_upspins_.push_back(state.spins);
      my_upspins_offset_[state.spins] = offset;
      offset += hs_holes_in_ups_.size();
    }
  }

  // Check whether every upspin configuration belongs to some process
  idx_t dim_test = my_upspins_.size();
  MPI_Allreduce(MPI_IN_PLACE, &dim_test, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  assert(dim_test == hs_upspins_.size());

  // Set dimensions
  local_dim_ = offset;
  MPI_Allreduce(&local_dim_, &dim_, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&local_dim_, &max_local_dim_, 1, MPI_UNSIGNED_LONG, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(&local_dim_, &min_local_dim_, 1, MPI_UNSIGNED_LONG, MPI_MIN,
                MPI_COMM_WORLD);

  // Collect all downspin configurations belonging to this mpi_rank
  offset = 0;
  for (auto state : hs_downspins_) {
    if (mpi_rank_of_spins(state.spins) == mpi_rank_) {
      my_downspins_.push_back(state.spins);
      my_downspins_offset_[state.spins] = offset;
      offset += hs_holes_in_downs_.size();
    }
  }
  local_dim_downspins_ = offset;

  // Check whether also local downspin dimensions sum up to the total dim
  idx_t dim_downspins;
  MPI_Allreduce(&local_dim_downspins_, &dim_downspins, 1, MPI_UNSIGNED_LONG,
                MPI_SUM, MPI_COMM_WORLD);
  assert(dim_downspins == dim_);

  ///////////////////////////////////////////////
  // Determine Communication patterns
  // for forward sending downspins
  n_downspins_i_send_forward_.resize(mpi_size_, 0);
  n_downspins_i_recv_forward_.resize(mpi_size_, 0);
  n_downspins_i_send_forward_offsets_.resize(mpi_size_, 0);
  n_downspins_i_recv_forward_offsets_.resize(mpi_size_, 0);
  for (auto holes : hs_holes_in_ups_) {
    for (auto upspins : my_upspins_) {
      auto downspins = up_hole_to_down(state_spinhalf<bit_t>{upspins}, holes);
      int destination_mpi_rank = mpi_rank_of_spins(downspins.spins);
      ++n_downspins_i_send_forward_[destination_mpi_rank];
    }
    // if (mpi_rank_ == 0) std::cout << std::endl;
  }

  // Communicate number of send/receive states
  lila::MPI_Alltoall<int>(n_downspins_i_send_forward_.data(), 1,
                          n_downspins_i_recv_forward_.data(), 1,
                          MPI_COMM_WORLD);

  // Compute offsets and total number of send/receive states
  for (int m = 0; m < mpi_size_; ++m)
    for (int n = 0; n < m; ++n) {
      n_downspins_i_send_forward_offsets_[m] += n_downspins_i_send_forward_[n];
      n_downspins_i_recv_forward_offsets_[m] += n_downspins_i_recv_forward_[n];
    }

  sum_n_downspins_i_send_forward_ =
      std::accumulate(n_downspins_i_send_forward_.begin(),
                      n_downspins_i_send_forward_.end(), 0);
  sum_n_downspins_i_recv_forward_ =
      std::accumulate(n_downspins_i_recv_forward_.begin(),
                      n_downspins_i_recv_forward_.end(), 0);

  // Communicate which downstates are sent/received by whom
  // std::vector<bit_t>
  // downspins_i_send_forward_(sum_n_downspins_i_send_forward_, 0);
  // std::vector<bit_t>
  // upspins_i_send_forward_(sum_n_downspins_i_send_forward_, 0);
  downspins_i_send_forward_ =
      std::vector<bit_t>(sum_n_downspins_i_send_forward_, 0);
  upspins_i_send_forward_ =
      std::vector<bit_t>(sum_n_downspins_i_send_forward_, 0);

  std::vector<idx_t> n_downspins_already_prepared(mpi_size_, 0);
  for (auto holes : hs_holes_in_ups_)
    for (auto upspins : my_upspins_) {
      auto downspins = up_hole_to_down(state_spinhalf<bit_t>{upspins}, holes);
      int destination_mpi_rank = mpi_rank_of_spins(downspins.spins);

      idx_t idx = n_downspins_i_send_forward_offsets_[destination_mpi_rank] +
                  n_downspins_already_prepared[destination_mpi_rank]++;
      downspins_i_send_forward_[idx] = downspins.spins;
      upspins_i_send_forward_[idx] = upspins;
    }
  downspins_i_recv_forward_.resize(sum_n_downspins_i_recv_forward_);
  upspins_i_recv_forward_.resize(sum_n_downspins_i_recv_forward_);

  lila::MPI_Alltoallv<bit_t>(
      downspins_i_send_forward_.data(), n_downspins_i_send_forward_.data(),
      n_downspins_i_send_forward_offsets_.data(),
      downspins_i_recv_forward_.data(), n_downspins_i_recv_forward_.data(),
      n_downspins_i_recv_forward_offsets_.data(), MPI_COMM_WORLD);
  lila::MPI_Alltoallv<bit_t>(
      upspins_i_send_forward_.data(), n_downspins_i_send_forward_.data(),
      n_downspins_i_send_forward_offsets_.data(),
      upspins_i_recv_forward_.data(), n_downspins_i_recv_forward_.data(),
      n_downspins_i_recv_forward_offsets_.data(), MPI_COMM_WORLD);

  // downspins_i_send_forward_ and upspins_i_send_forward_ not needed anymore
  // downspins_i_send_forward_.clear();
  // upspins_i_send_forward_.clear();
  // downspins_i_send_forward_.shrink_to_fit();
  // upspins_i_send_forward_.shrink_to_fit();

  ///////////////////////////////////////////////
  // Determine Communication patterns
  // for back sending downspins
  n_upspins_i_send_back_.resize(mpi_size_, 0);
  n_upspins_i_recv_back_.resize(mpi_size_, 0);
  n_upspins_i_send_back_offsets_.resize(mpi_size_, 0);
  n_upspins_i_recv_back_offsets_.resize(mpi_size_, 0);

  for (auto downspins : my_downspins_)
    for (auto holes : hs_holes_in_downs_) {
      bit_t upspins = down_hole_to_up(downspins, (bit_t)holes);
      int destination_mpi_rank = mpi_rank_of_spins(upspins);
      ++n_upspins_i_send_back_[destination_mpi_rank];
    }

  // Communicate number of send/receive states
  lila::MPI_Alltoall<int>(n_upspins_i_send_back_.data(), 1,
                          n_upspins_i_recv_back_.data(), 1, MPI_COMM_WORLD);

  // Compute offsets and total number of send/receive states
  for (int m = 0; m < mpi_size_; ++m)
    for (int n = 0; n < m; ++n) {
      n_upspins_i_send_back_offsets_[m] += n_upspins_i_send_back_[n];
      n_upspins_i_recv_back_offsets_[m] += n_upspins_i_recv_back_[n];
    }

  sum_n_upspins_i_send_back_ = std::accumulate(n_upspins_i_send_back_.begin(),
                                               n_upspins_i_send_back_.end(), 0);
  sum_n_upspins_i_recv_back_ = std::accumulate(n_upspins_i_recv_back_.begin(),
                                               n_upspins_i_recv_back_.end(), 0);

  // Communicate which upstates are sent/received by whom
  // std::vector<bit_t> downspins_i_send_back_(sum_n_upspins_i_send_back_, 0);
  // std::vector<bit_t> upspins_i_send_back_(sum_n_upspins_i_send_back_, 0);
  downspins_i_send_back_ = std::vector<bit_t>(sum_n_upspins_i_send_back_, 0);
  upspins_i_send_back_ = std::vector<bit_t>(sum_n_upspins_i_send_back_, 0);

  std::vector<idx_t> n_upspins_already_prepared(mpi_size_, 0);
  // if (mpi_rank_ == 0) std::cout << std::endl;

  for (auto holes : hs_holes_in_downs_) {
    for (auto downspins : my_downspins_) {
      bit_t upspins = down_hole_to_up(downspins, (bit_t)holes);
      int destination_mpi_rank = mpi_rank_of_spins(upspins);

      idx_t idx = n_upspins_i_send_back_offsets_[destination_mpi_rank] +
                  n_upspins_already_prepared[destination_mpi_rank]++;
      downspins_i_send_back_[idx] = downspins;
      upspins_i_send_back_[idx] = upspins;
    }
    // if (mpi_rank_ == 0) std::cout << std::endl;
  }

  downspins_i_recv_back_.resize(sum_n_upspins_i_recv_back_);
  upspins_i_recv_back_.resize(sum_n_upspins_i_recv_back_);

  lila::MPI_Alltoallv(
      downspins_i_send_back_.data(), n_upspins_i_send_back_.data(),
      n_upspins_i_send_back_offsets_.data(), downspins_i_recv_back_.data(),
      n_upspins_i_recv_back_.data(), n_upspins_i_recv_back_offsets_.data(),
      MPI_COMM_WORLD);
  lila::MPI_Alltoallv(
      upspins_i_send_back_.data(), n_upspins_i_send_back_.data(),
      n_upspins_i_send_back_offsets_.data(), upspins_i_recv_back_.data(),
      n_upspins_i_recv_back_.data(), n_upspins_i_recv_back_offsets_.data(),
      MPI_COMM_WORLD);

  // downspins_i_send_back_ and upspins_i_send_back_ not needed anymore
  // downspins_i_send_back_.clear();
  // upspins_i_send_back_.clear();
  // downspins_i_send_back_.shrink_to_fit();
  // upspins_i_send_back_.shrink_to_fit();

  // Allocate send/receive buffers
  buffer_size_ = std::max(
      std::max(sum_n_downspins_i_send_forward_,
               sum_n_downspins_i_recv_forward_),
      std::max(sum_n_upspins_i_send_back_, sum_n_upspins_i_recv_back_));

  try {
    send_buffer_.resize(buffer_size_);
    recv_buffer_.resize(buffer_size_);
  } catch (...) {
    std::cerr << "[ " << mpi_rank_
              << " ] Error: Could not allocate send/receive buffers!"
              << std::endl
              << std::flush;
    MPI_Abort(MPI_COMM_WORLD, 3);
  }

  // Create downspins_table_ and upspins_table
  for (auto upspins : my_upspins_)
    for (auto holes : hs_holes_in_ups_) {
      bit_t downspins = up_hole_to_down(upspins, (bit_t)holes);
      downspins_table_.push_back(downspins);
    }

  for (auto downspins : my_downspins_)
    for (auto holes : hs_holes_in_downs_) {
      bit_t upspins = down_hole_to_up(downspins, (bit_t)holes);
      upspins_table_.push_back(upspins);
    }
}

template class TJModelMPI<double, uint32>;
template class TJModelMPI<complex, uint32>;
template class TJModelMPI<double, uint64>;
template class TJModelMPI<complex, uint64>;

} // namespace hydra
