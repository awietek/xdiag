#include "hubbardmodelmpi.h"
#include "hubbardmodeldetail.h"
#include <hydra/utils/bitops.h>
#include <hydra/utils/complex.h>
#include <iostream>

namespace hydra {

template <class coeff_t, class bit_t, class idx_t>
HubbardModelMPI<coeff_t, bit_t, idx_t>::HubbardModelMPI(BondList bondlist,
                                                        Couplings couplings,
                                                        qn_electron qn)
    : n_sites_(bondlist.n_sites()), qn_(qn) {
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);


  detail::set_hubbard_terms<coeff_t>(
      bondlist, couplings, hoppings_, hopping_amplitudes_, currents_,
      current_amplitudes_, interactions_, interaction_strengths_, onsites_,
      onsite_potentials_, szszs_, szsz_amplitudes_, exchanges_,
      exchange_amplitudes_, U_);

  initialize();
}

template <class coeff_t, class bit_t, class idx_t>
void HubbardModelMPI<coeff_t, bit_t, idx_t>::set_qn(qn_electron qn) {
  qn_ = qn;
  initialize();
}

template <class coeff_t, class bit_t, class idx_t>
void HubbardModelMPI<coeff_t, bit_t, idx_t>::apply_hamiltonian(
    const lila::VectorMPI<coeff_t> &in_vec, lila::VectorMPI<coeff_t> &out_vec,
    bool verbose) {
  using utils::gbit;
  using utils::gbits;
  using utils::popcnt;

  if ((mpi_rank_ == 0) && verbose) {
    printf("max_local_dim: %lu\n", max_local_dim_);
    printf("min_local_dim: %lu\n", min_local_dim_);
  }
  Zeros(out_vec.vector_local());

  double t1 = MPI_Wtime();

  // Apply Diagonal U terms
  if (std::abs(U_) > 1e-14) {
    idx_t upspin_idx = 0;
    for (const bit_t &upspins : my_upspins_) // loop over upspins of process
    {
      idx_t upspin_offset = my_upspins_offset_[upspins];
      idx_t downspin_offset = 0;
      for (auto downspins : basis_dn_) // loop over all downspins
      {
        idx_t idx = upspin_offset + downspin_offset;
        double coeff = U_ * (double)popcnt(upspins & downspins.spins);
        out_vec.vector_local()(idx) += coeff * in_vec.vector_local()(idx);
        ++downspin_offset;
      }
      ++upspin_idx;
    }
  }

  // Apply Diagonal V terms
  int interaction_idx = 0;
  for (auto pair : interactions_) {
    const int s1 = pair.first;
    const int s2 = pair.second;
    const double V = interaction_strengths_[interaction_idx];

    if (std::abs(V) > 1e-14) {
      idx_t upspin_idx = 0;
      for (const bit_t &upspins : my_upspins_) // loop over upspins of process
      {
        idx_t upspin_offset = my_upspins_offset_[upspins];
        idx_t downspin_offset = 0;
        for (auto downspins : basis_dn_) // loop over all downspins
        {
          idx_t idx = upspin_offset + downspin_offset;
          auto coeff = V * (double)((gbit(upspins, s1) + gbit(downspins.spins, s1)) *
                                    (gbit(upspins, s2) + gbit(downspins.spins, s2)));
          out_vec.vector_local()(idx) += coeff * in_vec(idx);
          ++downspin_offset;
        }
        ++upspin_idx;
      }
    } // (std::abs(V) > 1e-14)

    ++interaction_idx;
  }

  // Apply onsite chemical potential
  int onsite_idx = 0;
  for (auto site : onsites_) {
    const double mu = onsite_potentials_[onsite_idx];

    if (std::abs(mu) > 1e-14) {
      idx_t upspin_idx = 0;
      for (const bit_t &upspins : my_upspins_) // loop over upspins of process
      {
        idx_t upspin_offset = my_upspins_offset_[upspins];
        idx_t downspin_offset = 0;
        for (auto downspins : basis_dn_) // loop over all downspins
        {
          idx_t idx = upspin_offset + downspin_offset;
          auto coeff =
              mu * (double)((gbit(upspins, site) + gbit(downspins.spins, site)));
          out_vec.vector_local()(idx) -= coeff * in_vec(idx);
          ++downspin_offset;
        }
        ++upspin_idx;
      }
    } // (std::abs(mu) > 1e-14)

    ++onsite_idx;
  }

  // Apply szsz interactions
  int szsz_idx = 0;
  for (auto pair : szszs_) {
    const int s1 = std::min(pair.first, pair.second);
    const int s2 = std::max(pair.first, pair.second);
    const double jz = szsz_amplitudes_[szsz_idx] * 0.25;
    if (std::abs(jz) > 1e-14) {
      idx_t upspin_idx = 0;
      for (const bit_t &upspins : my_upspins_) {
        idx_t upspin_offset = my_upspins_offset_[upspins];
        idx_t downspin_offset = 0;
        for (auto downspins : basis_dn_) {
          idx_t idx = upspin_offset + downspin_offset;
          auto coeff =
              jz * (((double)gbit(upspins, s1) - (double)gbit(downspins.spins, s1)) *
                    ((double)gbit(upspins, s2) - (double)gbit(downspins.spins, s2)));
          out_vec(idx) += coeff * in_vec(idx);
          ++downspin_offset;
        }
        ++upspin_idx;
      }
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

    // Find out how many states is send to each process
    std::vector<int> n_states_i_send(mpi_size_, 0);
    if (std::abs(jx) > 1e-14) {
      // Flip states and check out how much needs to be communicated
      for (const bit_t &upspins : my_upspins_)
        for (auto downspins : basis_dn_) {
          if ((popcnt(upspins & flipmask) == 1) &&
              (popcnt(downspins.spins & flipmask) == 1) &&
              popcnt((downspins.spins & flipmask) & (upspins & flipmask)) == 0) {
            bit_t flipped_upspins = upspins ^ flipmask;
            int target = mpi_rank_of_spins(flipped_upspins);
            ++n_states_i_send[target];
          }
        }

      // printf("n_states_i_send\n");
      // for (int m=0; m<mpi_size_; ++m)
      // 	printf(" [%d] -> [%d] %d\n", mpi_rank_, m, n_states_i_send[m]);

      // Exchange information on who sends how much to whom
      std::vector<int> n_states_i_recv(mpi_size_, 0);
      lila::MPI_Alltoall<int>(n_states_i_send.data(), 1, n_states_i_recv.data(),
                              1, MPI_COMM_WORLD);

      // printf("n_states_i_recv\n");
      // for (int m=0; m<mpi_size_; ++m)
      // 	printf(" [%d] <- [%d] %d\n", mpi_rank_, m, n_states_i_recv[m]);

      // Sum up states sent/recvd and eventually resize buffers
      idx_t sum_n_states_i_send = std::accumulate(
          n_states_i_send.begin(), n_states_i_send.end(), (idx_t)0);
      idx_t sum_n_states_i_recv = std::accumulate(
          n_states_i_recv.begin(), n_states_i_recv.end(), (idx_t)0);

      if (sum_n_states_i_send > send_buffer_.size())
        send_buffer_.resize(sum_n_states_i_send);
      if (sum_n_states_i_recv > recv_buffer_.size())
        recv_buffer_.resize(sum_n_states_i_recv);

      // printf("sum_n_states_i_send: %ld\n", sum_n_states_i_send);
      // printf("sum_n_states_i_recv: %ld\n", sum_n_states_i_recv);

      // Compute offsets of send/receive states
      std::vector<int> n_states_i_send_offsets(mpi_size_, 0);
      std::vector<int> n_states_i_recv_offsets(mpi_size_, 0);
      for (int m = 0; m < mpi_size_; ++m)
        for (int n = 0; n < m; ++n) {
          n_states_i_send_offsets[m] += n_states_i_send[n];
          n_states_i_recv_offsets[m] += n_states_i_recv[n];
        }

      // printf("n_states_i_send_offsets\n");
      // for (int m=0; m<mpi_size_; ++m)
      // 	printf(" [%d] -> [%d] %d\n", mpi_rank_, m,
      // n_states_i_send_offsets[m]); printf("n_states_i_recv_offsets\n"); for
      // (int m=0; m<mpi_size_; ++m) 	printf(" [%d] <- [%d] %d\n", mpi_rank_,
      // m, n_states_i_recv_offsets[m]);

      // Flip states and check out how much needs to be communicated
      std::vector<int> n_states_prepared(mpi_size_, 0);
      for (const bit_t &upspins : my_upspins_) {
        idx_t upspin_offset = my_upspins_offset_[upspins];
        idx_t downspin_offset = 0;
        for (auto downspins : basis_dn_) {

          if ((popcnt(upspins & flipmask) == 1) &&
              (popcnt(downspins.spins & flipmask) == 1) &&
              popcnt((downspins.spins & flipmask) & (upspins & flipmask)) == 0) {
            idx_t idx = upspin_offset + downspin_offset;
            bit_t flipped_upspins = upspins ^ flipmask;
            int target = mpi_rank_of_spins(flipped_upspins);
            int send_idx =
                n_states_i_send_offsets[target] + n_states_prepared[target];

            send_buffer_[send_idx] = in_vec.vector_local()(idx);
            ++n_states_prepared[target];
          }
          ++downspin_offset;
        }
      }

      // LilaPrint(in_vec.vector_local());
      // for (int i=0; i<sum_n_states_i_send; ++i)
      // 	printf("send_buf[%d] = %f\n", i, send_buffer_[i]);

      // Check, whether correct number of states has been prepared
      for (int m = 0; m < mpi_size_; ++m)
        assert(n_states_prepared[m] == n_states_i_send[m]);

      // Alltoall call
      lila::MPI_Alltoallv<coeff_t>(send_buffer_.data(), n_states_i_send.data(),
                                   n_states_i_send_offsets.data(),
                                   recv_buffer_.data(), n_states_i_recv.data(),
                                   n_states_i_recv_offsets.data(),
                                   MPI_COMM_WORLD);

      // for (int i=0; i<sum_n_states_i_recv; ++i)
      // 	printf("recv_buf[%d] = %f\n", i, recv_buffer_[i]);

      // Get the original upspin configuration and its source proc
      std::vector<std::vector<bit_t>> upspins_i_get_from_proc(mpi_size_);
      for (const bit_t &upspins : my_upspins_) {
        bit_t flipped_upspins = upspins ^ flipmask;
        int source = mpi_rank_of_spins(flipped_upspins);
        upspins_i_get_from_proc[source].push_back(upspins);
      }

      // for (int i=0; i<mpi_size_; ++i)
      // 	for (auto upspins : upspins_i_get_from_proc[i])
      // 	  printf("upspins_i_get: %d\n", upspins);

      // Sort the upspins and add coefficients to outvec
      // LilaPrint(out_vec.vector_local());

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
          idx_t upspin_offset = my_upspins_offset_[upspins];
          idx_t downspin_offset = 0;

          // Skip upspins if exchange isn't possible
          if (popcnt(upspins & flipmask) != 1)
            continue;

          bit_t downmask =
              gbit(upspins, s1) ? (bit_t)1 << s2 : (bit_t)1 << s1;

          for (auto downspins : basis_dn_) {
            if ((downspins.spins & flipmask) == downmask) {
              idx_t target_idx = upspin_offset + downspin_offset;
              out_vec.vector_local()(target_idx) += jx * recv_buffer_[recv_idx];
              ++recv_idx;
            }
            ++downspin_offset;
          }
        }
      }
      // LilaPrint(out_vec.vector_local());
    }
    ++exchange_idx;
  }
  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  exch: %3.4f\n", t2 - t1);

  // Apply hoppings on downspins
  double ht1 = MPI_Wtime();
  t1 = MPI_Wtime();
  int hopping_idx = 0;
  for (auto pair : hoppings_) {
    const int s1 = std::min(pair.first, pair.second);
    const int s2 = std::max(pair.first, pair.second);
    const coeff_t t = hopping_amplitudes_[hopping_idx];
    const bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    const bit_t secondmask = (bit_t)1 << s2;

    if (std::abs(t) > 1e-14) {
      // Loop over all configurations
      idx_t upspin_idx = 0;
      for (const bit_t &upspins : my_upspins_) {
        idx_t upspin_offset = my_upspins_offset_[upspins];
        idx_t downspin_offset = 0;
        for (auto downspins : basis_dn_) {
          // Check if hopping is possible
          if (((downspins.spins & flipmask) != 0) &&
              ((downspins.spins & flipmask) != flipmask)) {
            const double fermi =
                popcnt(gbits(upspins ^ (downspins.spins ^ secondmask), s2 - s1,
                             s1 + 1)) &
                        1
                    ? 1.
                    : -1.;

            idx_t idx = upspin_offset + downspin_offset;
            bit_t new_downspins = downspins.spins ^ flipmask;
            idx_t new_idx = upspin_offset + index_dn_.index({new_downspins});

            if (downspins.spins & secondmask) {
              out_vec.vector_local()(new_idx) +=
                  fermi * t * in_vec.vector_local()(idx);
            } else {
              out_vec.vector_local()(new_idx) -=
                  fermi * lila::conj(t) * in_vec.vector_local()(idx);
            }
          }

          ++downspin_offset;
        }
        ++upspin_idx;
      }
    }
    ++hopping_idx;
  } // hopping on downspin

  // Apply currents on downspins
  int current_idx = 0;
  for (auto pair : currents_) {
    const int s1 = std::min(pair.first, pair.second);
    const int s2 = std::max(pair.first, pair.second);
    const coeff_t t = current_amplitudes_[current_idx];
    const bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    // printf("down p1: %d, p2: %d, s1: %d, s2: %d\n",
    // 	 pair.first, pair.second, s1, s2);

    if (std::abs(t) > 1e-14) {
      // Loop over all configurations
      idx_t upspin_idx = 0;
      for (const bit_t &upspins : my_upspins_) {
        idx_t upspin_offset = my_upspins_offset_[upspins];
        idx_t downspin_offset = 0;
        for (auto downspins : basis_dn_) {

          // Check if current is possible
          if (((downspins.spins & flipmask) != 0) &&
              ((downspins.spins & flipmask) != flipmask)) {
            double fermi =
                popcnt(gbits(downspins.spins, s2 - s1 - 1, s1 + 1)) % 2 == 0 ? 1.
                                                                       : -1.;
            double dir =
                gbit(downspins.spins, pair.first) ? 1. : -1.; // use pair, not s1 s2

            idx_t idx = upspin_offset + downspin_offset;
            bit_t new_downspins = downspins.spins ^ flipmask;
            idx_t new_idx = upspin_offset + index_dn_.index({new_downspins});
            out_vec.vector_local()(new_idx) -=
                fermi * dir * t * in_vec.vector_local()(idx);
          }

          ++downspin_offset;
        }
        ++upspin_idx;
      }
    }
    ++hopping_idx;
  } // hopping on downspin

  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  down: %3.4f\n", t2 - t1);

  //
  // Communication for hoppings on upspins
  //

  // Send configurations
  t1 = MPI_Wtime();
  std::vector<idx_t> n_states_already_prepared(mpi_size_, 0);
  idx_t downspin_offset = 0;
  for (auto downspins : basis_dn_) {
    int destination_mpi_rank = mpi_rank_of_spins(downspins.spins);
    for (idx_t upspin_idx = 0; upspin_idx < my_upspins_.size(); ++upspin_idx) {
      idx_t send_idx =
          n_downspins_i_send_forward_offsets_[destination_mpi_rank] +
          n_states_already_prepared[destination_mpi_rank]++;
      idx_t upspin_offset = my_upspins_offset_[my_upspins_[upspin_idx]];
      idx_t idx = upspin_offset + downspin_offset;
      send_buffer_[send_idx] = in_vec.vector_local()(idx);
    }
    ++downspin_offset;
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
  // 	for (const bit_t& upspins : my_upspins_)  // loop over upspins of
  // process
  // 	  {
  // 	    idx_t upspin_offset = my_upspins_offset_[upspins];
  // 	    idx_t downspin_offset=0;
  // 	    for (bit_t downspins : basis_dn_)  // loop over all downspins
  // 	      {
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
  // 	printf("[%d]:\n", mpi_rank_);

  // 	for (int sendto = 0; sendto < mpi_size_; ++sendto)
  // 	  {
  // 	    int from = n_downspins_i_send_forward_offsets_[sendto];
  // 	    int to = (sendto == mpi_size_-1) ? sum_n_downspins_i_send_forward_ :
  // n_downspins_i_send_forward_offsets_[sendto+1]; 	    for (int idx = from;
  // idx < to; ++idx) 	      printf("%f -> [%d] %s; %s\n",
  // send_buffer_[idx], sendto, 		     PrintSpinhalf(n_sites_,
  // upspins_i_send_forward_[idx]).c_str(), PrintSpinhalf(n_sites_,
  // downspins_i_send_forward_[idx]).c_str());
  // 	  }
  // 	printf("\n");
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
  for (idx_t idx = 0; idx < sum_n_downspins_i_recv_forward_; ++idx) {
    bit_t upspins = upspins_i_recv_forward_[idx];
    bit_t downspins = downspins_i_recv_forward_[idx];

    idx_t sorted_idx =
      my_downspins_offset_[downspins] + index_up_.index({upspins});
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
  // 	for (const bit_t& downspins : my_downspins_)
  // 	  {
  // 	    for (bit_t upspins : basis_up_)
  // 	      {
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

  // Apply upspin hoppings on the sorted send buffer
  // write the results to the receive buffer

  // // DEBUG START
  // // JUST COPY INVEC FOR DEBUGGING
  // std::copy(send_buffer_.begin(), send_buffer_.end(), recv_buffer_.begin());
  // // DEBUG END

  std::fill(recv_buffer_.begin(), recv_buffer_.end(), 0);

  t1 = MPI_Wtime();

  // Apply hoppings on upspins
  hopping_idx = 0;
  for (auto pair : hoppings_) {
    const int s1 = std::min(pair.first, pair.second);
    const int s2 = std::max(pair.first, pair.second);
    const coeff_t t = hopping_amplitudes_[hopping_idx];
    const bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    const bit_t secondmask = (bit_t)1 << s2;

    if (std::abs(t) > 1e-14) {
      // Loop over all configurations
      idx_t downspin_idx = 0;
      for (const bit_t &downspins : my_downspins_) {
        idx_t downspin_offset = my_downspins_offset_[downspins];
        idx_t upspin_offset = 0;
        for (auto upspins : basis_up_) {
          // Check if hopping is possible
          if (((upspins.spins & flipmask) != 0) &&
              ((upspins.spins & flipmask) != flipmask)) {
            const double fermi =
                popcnt(gbits(upspins.spins ^ downspins, s2 - s1, s1)) & 1 ? 1. : -1.;

            idx_t idx = upspin_offset + downspin_offset;
            bit_t new_upspins = upspins.spins ^ flipmask;
            idx_t new_idx = downspin_offset + index_up_.index({new_upspins});
            if (upspins.spins & secondmask) {
              recv_buffer_[new_idx] += fermi * t * send_buffer_[idx];
            } else {
              recv_buffer_[new_idx] -=
                  fermi * lila::conj(t) * send_buffer_[idx];
            }
          }
          ++upspin_offset;
        }
        ++downspin_idx;
      }
    }
    ++hopping_idx;
  } // hopping on upspins

  // Apply currents on upspins
  current_idx = 0;
  for (auto pair : currents_) {
    const int s1 = std::min(pair.first, pair.second);
    const int s2 = std::max(pair.first, pair.second);
    const coeff_t t = current_amplitudes_[current_idx];
    const bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    // printf("up p1: %d, p2: %d, s1: %d, s2: %d\n",
    // 	 pair.first, pair.second, s1, s2);
    if (std::abs(t) > 1e-14) {
      // Loop over all configurations
      idx_t downspin_idx = 0;
      for (const bit_t &downspins : my_downspins_) {
        idx_t downspin_offset = my_downspins_offset_[downspins];
        idx_t upspin_offset = 0;
        for (auto upspins : basis_up_) {
          // Check if current is possible
          if (((upspins.spins & flipmask) != 0) &&
              ((upspins.spins & flipmask) != flipmask)) {
            double fermi =
                popcnt(gbits(upspins.spins, s2 - s1 - 1, s1 + 1)) % 2 == 0 ? 1. : -1.;
            double dir =
                gbit(upspins.spins, pair.first) ? 1. : -1.; // use pair, not s1 s2
            // printf("up p1: %d, p2: %d, s1: %d, s2: %d, v1: %d, v2: %d, dir:
            // %f\n", 	 pair.first, pair.second, s1, s2, gbit(upspins,
            // pair.first), gbit(upspins, pair.second), dir);

            idx_t idx = upspin_offset + downspin_offset;
            bit_t new_upspins = upspins.spins ^ flipmask;
            idx_t new_idx = downspin_offset + index_up_.index({new_upspins});
            recv_buffer_[new_idx] -= fermi * dir * t * send_buffer_[idx];
          }

          ++upspin_offset;
        }
        ++downspin_idx;
      }
    }
    ++current_idx;
  } // current on upspins

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
  // 	for (const bit_t& downspins : my_downspins_)
  // 	  {
  // 	    for (bit_t upspins : basis_up_)
  // 	      {
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
  idx_t upspin_offset = 0;
  for (auto upspins : basis_up_) {
    int destination_mpi_rank = mpi_rank_of_spins(upspins.spins);
    for (idx_t downspin_idx = 0; downspin_idx < my_downspins_.size();
         ++downspin_idx) {
      idx_t send_idx = n_upspins_i_send_back_offsets_[destination_mpi_rank] +
                        n_states_already_prepared[destination_mpi_rank]++;
      idx_t downspin_offset =
          my_downspins_offset_[my_downspins_[downspin_idx]];
      idx_t idx = upspin_offset + downspin_offset;
      send_buffer_[send_idx] = recv_buffer_[idx];
    }
    ++upspin_offset;
  }
  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  send back (prepare):   %3.4f\n", t2 - t1);

  // // DEBUG START
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(mpi_rank_ == 0) printf("BACK SEND BUFFER
  // --------------------------------\n"); MPI_Barrier(MPI_COMM_WORLD); for(int
  // nt = 0; nt < mpi_size_; ++nt)
  //   {
  //     if(mpi_rank_ == nt)
  //       {
  // 	printf("[%d]:\n", mpi_rank_);

  // 	for (int sendto = 0; sendto < mpi_size_; ++sendto)
  // 	  {
  // 	    int from = n_upspins_i_send_back_offsets_[sendto];
  // 	    int to = (sendto == mpi_size_-1) ? sum_n_upspins_i_send_back_ :
  // n_upspins_i_send_back_offsets_[sendto+1]; 	    for (int idx = from; idx <
  // to;
  // ++idx) 	      printf("%f -> [%d] %s; %s \n", send_buffer_[idx], sendto,
  // 		     PrintSpinhalf(n_sites_, upspins_i_send_back_[idx]).c_str(),
  // 		     PrintSpinhalf(n_sites_,
  // downspins_i_send_back_[idx]).c_str());
  // 	  }
  // 	printf("\n");
  //       }
  //     MPI_Barrier(MPI_COMM_WORLD);
  //   }
  // MPI_Barrier(MPI_COMM_WORLD);
  // // DEBUG END

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
  for (idx_t idx = 0; idx < sum_n_upspins_i_recv_back_; ++idx) {
    bit_t upspins = upspins_i_recv_back_[idx];
    bit_t downspins = downspins_i_recv_back_[idx];
    idx_t sorted_idx =
      my_upspins_offset_[upspins] + index_dn_.index({downspins});
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
  // 	for (const bit_t& upspins : my_upspins_)
  // 	  {
  // 	    for (bit_t downspins : basis_dn_)
  // 	      {
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
  // printf("[%d] send_buffer_size: %d, in_vec.vector_local().size(): %d,
  // out_vec.vector_local().size(): %d\n", 	     mpi_rank_,
  // send_buffer_.size(), in_vec.vector_local().size(),
  // out_vec.vector_local().size());

  for (int k = 0; k < in_vec.vector_local().size(); ++k)
    out_vec.vector_local()(k) += send_buffer_[k];
  t2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  fill outvec:   %3.4f\n", t2 - t1);

  double ht2 = MPI_Wtime();
  if ((mpi_rank_ == 0) && verbose)
    printf("  hopp: %3.4f\n", ht2 - ht1);
}

template <class coeff_t, class bit_t, class idx_t>
qn_electron HubbardModelMPI<coeff_t, bit_t, idx_t>::apply_fermion(
    const lila::VectorMPI<coeff_t> &in_vec, lila::VectorMPI<coeff_t> &out_vec,
    std::string type, int site) {
  using utils::gbit;
  using utils::gbits;
  using utils::popcnt;
  assert(site < n_sites_);

  auto qn_after = qn_;
  if (type == "cdagup") {
    ++qn_after.n_up;
    if (mpi_rank_ == 0)
      printf("cdagup not implemented! Use cdagdn instead.\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else if (type == "cup") {
    --qn_after.n_up;
    if (mpi_rank_ == 0)
      printf("cup not implemented! Use cdn instead.\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else if (type == "cdagdn")
    ++qn_after.n_dn;
  else if (type == "cdn")
    --qn_after.n_dn;
  else {
    if (mpi_rank_ == 0)
      printf("Error in apply_fermion: Invalid fermion type!\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  BasisSpinHalf<bit_t> basis_dn_after(n_sites_, qn_after.n_dn);
  IndexTable<BasisSpinHalf<bit_t>, idx_t> index_dn_after(basis_dn_after);

  // Compute new offsets for altered hilbertspace
  std::unordered_map<bit_t, idx_t> my_upspins_offset_after;
  for (auto old_offset : my_upspins_offset_) {
    bit_t state = old_offset.first;
    idx_t off = old_offset.second;
    idx_t new_off = off / index_dn_.size() * index_dn_after.size();
    my_upspins_offset_after[state] = new_off;
  }

  // Allocate out_vec
  idx_t local_dim_after = my_upspins_.size() * index_dn_after.size();
  try {
    out_vec.resize(local_dim_after);
  } catch (...) {
    std::cerr << "[ " << mpi_rank_
              << " ] Error: Could not allocate fermion vector!" << std::endl
              << std::flush;
    MPI_Abort(MPI_COMM_WORLD, 4);
  }

  // Loop over all configurations
  idx_t upspin_idx = 0;
  for (const bit_t &upspins : my_upspins_) {
    idx_t upspin_offset = my_upspins_offset_[upspins];
    idx_t upspin_offset_after = my_upspins_offset_after[upspins];
    idx_t downspin_offset = 0;

    // Apply the electron part
    if (type == "cdagdn") {
      const bit_t sitemask = ((bit_t)1 << site);

      for (auto downspins : basis_dn_) {

        // raise local site val if 0
        if (gbit(downspins.spins, site) == 0) {
          idx_t idx = upspin_offset + downspin_offset;
          bit_t new_downspins = downspins.spins;
          new_downspins |= sitemask;
          idx_t new_idx =
	    upspin_offset_after + index_dn_after.index({new_downspins});
          double fermi = popcnt(gbits(downspins.spins, site, 0)) % 2 == 0 ? 1. : -1.;
          out_vec.vector_local()(new_idx) += fermi * in_vec.vector_local()(idx);
        }

        ++downspin_offset;
      }
    } // cdagdn

    // Apply the hole part
    else if (type == "cdn") {
      const bit_t antisitemask = ~((bit_t)1 << site);

      for (auto downspins : basis_dn_) {

        // lower local site val if 1
        if (gbit(downspins.spins, site) == 1) {
          idx_t idx = upspin_offset + downspin_offset;
          bit_t new_downspins = downspins.spins;
          new_downspins &= antisitemask;
          idx_t new_idx =
	    upspin_offset_after + index_dn_after.index({new_downspins});
          double fermi = popcnt(gbits(downspins.spins, site, 0)) % 2 == 0 ? 1. : -1.;
          out_vec.vector_local()(new_idx) += fermi * in_vec.vector_local()(idx);
        }

        ++downspin_offset;
      }
    } // cdn

    ++upspin_idx;
  } // loop over upspin configurations

  return qn_after;
}

template <class coeff_t, class bit_t, class idx_t>
void HubbardModelMPI<coeff_t, bit_t, idx_t>::initialize() {

  basis_up_ = BasisSpinHalf<bit_t>(n_sites_, qn_.n_up);
  basis_dn_ = BasisSpinHalf<bit_t>(n_sites_, qn_.n_dn);

  index_up_ = IndexTable<BasisSpinHalf<bit_t>, idx_t>(basis_up_);
  index_dn_ = IndexTable<BasisSpinHalf<bit_t>, idx_t>(basis_dn_);

  // clear previous members
  my_upspins_.clear();
  my_upspins_offset_.clear();
  my_downspins_.clear();
  my_downspins_offset_.clear();
  n_downspins_i_send_forward_.clear();
  n_downspins_i_recv_forward_.clear();
  n_downspins_i_send_forward_offsets_.clear();
  n_downspins_i_recv_forward_offsets_.clear();

  downspins_i_recv_forward_.clear();
  upspins_i_recv_forward_.clear();
  downspins_i_send_forward_offsets_.clear();
  downspins_i_recv_forward_offsets_.clear();

  sum_n_downspins_i_send_forward_ = 0;
  sum_n_downspins_i_recv_forward_ = 0;

  n_upspins_i_send_back_.clear();
  n_upspins_i_recv_back_.clear();
  n_upspins_i_send_back_offsets_.clear();
  n_upspins_i_recv_back_offsets_.clear();
  downspins_i_recv_back_.clear();
  upspins_i_recv_back_.clear();
  upspins_i_send_back_offsets_.clear();
  upspins_i_recv_back_offsets_.clear();

  sum_n_upspins_i_send_back_ = 0;
  sum_n_upspins_i_recv_back_ = 0;

  buffer_size_ = 0;
  send_buffer_.clear();
  recv_buffer_.clear();

  // Collect all upspin configurations belonging to this mpi_rank
  idx_t offset = 0;
  for (auto state : basis_up_) {
    if (mpi_rank_of_spins(state.spins) == mpi_rank_) {
      my_upspins_.push_back(state.spins);
      my_upspins_offset_[state.spins] = offset;
      offset += basis_dn_.size();
    }
  }

  // Check whether every upspin configuration belongs to some process
  idx_t dim_test = my_upspins_.size();
  MPI_Allreduce(MPI_IN_PLACE, &dim_test, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  assert(dim_test == basis_up_.size());

  local_dim_ = offset;
  MPI_Allreduce(&local_dim_, &dim_, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&local_dim_, &max_local_dim_, 1, MPI_UNSIGNED_LONG, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(&local_dim_, &min_local_dim_, 1, MPI_UNSIGNED_LONG, MPI_MIN,
                MPI_COMM_WORLD);

  // Collect all downspin configurations belonging to this mpi_rank
  offset = 0;
  for (auto state : basis_dn_) {
    if (mpi_rank_of_spins(state.spins) == mpi_rank_) {
      my_downspins_.push_back(state.spins);
      my_downspins_offset_[state.spins] = offset;
      offset += basis_up_.size();
    }
  }
  local_dim_downspins_ = offset;
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

  for (auto downspins : basis_dn_) {
    int destination_mpi_rank = mpi_rank_of_spins(downspins.spins);
    n_downspins_i_send_forward_[destination_mpi_rank] += my_upspins_.size();
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
  std::vector<bit_t> downspins_i_send_forward_(
      sum_n_downspins_i_send_forward_, 0);
  std::vector<bit_t> upspins_i_send_forward_(sum_n_downspins_i_send_forward_,
                                               0);

  // downspins_i_send_forward_ =
  // std::vector<bit_t>(sum_n_downspins_i_send_forward_, 0);
  // upspins_i_send_forward_=
  // std::vector<bit_t>(sum_n_downspins_i_send_forward_, 0);

  std::vector<idx_t> n_downspins_already_prepared(mpi_size_, 0);
  for (auto downspins : basis_dn_) {
    int destination_mpi_rank = mpi_rank_of_spins(downspins.spins);
    for (const bit_t &upspins : my_upspins_) // loop over upspins of process
    {
      idx_t idx = n_downspins_i_send_forward_offsets_[destination_mpi_rank] +
                   n_downspins_already_prepared[destination_mpi_rank]++;
      downspins_i_send_forward_[idx] = downspins.spins;
      upspins_i_send_forward_[idx] = upspins;
    }
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

  // downspins_i_send_forward_ and upspins_i_send_forward_ not neede anymore
  downspins_i_send_forward_.clear();
  upspins_i_send_forward_.clear();
  downspins_i_send_forward_.shrink_to_fit();
  upspins_i_send_forward_.shrink_to_fit();

  // // Print forward communication patterns
  // MPI_Barrier(MPI_COMM_WORLD);
  // for(int nt = 0; nt < mpi_size_; ++nt)
  //   {
  //     if(mpi_rank_ == nt)
  //       {
  // 	printf("[%d]: forward\n", mpi_rank_);
  // 	printf("my_upspins.size(): %d\n", my_upspins_.size());
  // 	printf("sum_n_downspins_i_send_forward_: %d\n",
  // sum_n_downspins_i_send_forward_); 	printf("sum_n_downspins_i_recv_forward_:
  // %d\n", sum_n_downspins_i_recv_forward_);

  // 	for(int k = 0; k < mpi_size_; ++k)
  // 	  {
  // 	    printf("-> %d (down): %d\n", k, n_downspins_i_send_forward_[k]);
  // 	    int from = n_downspins_i_send_forward_offsets_[k];
  // 	    int to = (k == mpi_size_-1) ? sum_n_downspins_i_send_forward_ :
  // n_downspins_i_send_forward_offsets_[k+1]; 	    for (int l = from; l < to;
  // ++l) 	      printf("  %s \n", PrintSpinhalf(n_sites_,
  // downspins_i_send_forward_[l]).c_str());
  // 	  }

  // 	for(int k = 0; k < mpi_size_; ++k)
  // 	  {
  // 	    printf("<- %d (down): %d\n", k, n_downspins_i_recv_forward_[k]);
  // 	    int from = n_downspins_i_recv_forward_offsets_[k];
  // 	    int to = (k == mpi_size_-1) ? sum_n_downspins_i_recv_forward_ :
  // n_downspins_i_recv_forward_offsets_[k+1]; 	    for (int l = from; l < to;
  // ++l) 	      printf("  %s \n", PrintSpinhalf(n_sites_,
  // downspins_i_recv_forward_[l]).c_str());
  // 	  }
  // 	printf("\n");
  //       }
  //     MPI_Barrier(MPI_COMM_WORLD);
  //   }
  // MPI_Barrier(MPI_COMM_WORLD);

  ///////////////////////////////////////////////
  // Determine Communication patterns
  // for back sending downspins
  n_upspins_i_send_back_.resize(mpi_size_, 0);
  n_upspins_i_recv_back_.resize(mpi_size_, 0);
  n_upspins_i_send_back_offsets_.resize(mpi_size_, 0);
  n_upspins_i_recv_back_offsets_.resize(mpi_size_, 0);

  for (auto upspins : basis_up_) {
    int destination_mpi_rank = mpi_rank_of_spins(upspins.spins);
    n_upspins_i_send_back_[destination_mpi_rank] += my_downspins_.size();
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
  std::vector<bit_t> downspins_i_send_back_(sum_n_upspins_i_send_back_, 0);
  std::vector<bit_t> upspins_i_send_back_(sum_n_upspins_i_send_back_, 0);

  // downspins_i_send_back_ = std::vector<bit_t>(sum_n_upspins_i_send_back_,
  // 0); upspins_i_send_back_ = std::vector<bit_t>(sum_n_upspins_i_send_back_,
  // 0);

  std::vector<idx_t> n_upspins_already_prepared(mpi_size_, 0);
  for (auto upspins : basis_up_) {
    int destination_mpi_rank = mpi_rank_of_spins(upspins.spins);
    for (const bit_t &downspins :
         my_downspins_) // loop over upspins of process
    {
      idx_t idx = n_upspins_i_send_back_offsets_[destination_mpi_rank] +
                   n_upspins_already_prepared[destination_mpi_rank]++;
      downspins_i_send_back_[idx] = downspins;
      upspins_i_send_back_[idx] = upspins.spins;
    }
  }
  downspins_i_recv_back_.resize(sum_n_upspins_i_recv_back_);
  upspins_i_recv_back_.resize(sum_n_upspins_i_recv_back_);

  lila::MPI_Alltoallv<bit_t>(
      downspins_i_send_back_.data(), n_upspins_i_send_back_.data(),
      n_upspins_i_send_back_offsets_.data(), downspins_i_recv_back_.data(),
      n_upspins_i_recv_back_.data(), n_upspins_i_recv_back_offsets_.data(),
      MPI_COMM_WORLD);
  lila::MPI_Alltoallv<bit_t>(
      upspins_i_send_back_.data(), n_upspins_i_send_back_.data(),
      n_upspins_i_send_back_offsets_.data(), upspins_i_recv_back_.data(),
      n_upspins_i_recv_back_.data(), n_upspins_i_recv_back_offsets_.data(),
      MPI_COMM_WORLD);

  // downspins_i_send_back_ and upspins_i_send_back_ not neede anymore
  downspins_i_send_back_.clear();
  upspins_i_send_back_.clear();
  downspins_i_send_back_.shrink_to_fit();
  upspins_i_send_back_.shrink_to_fit();

  // // Print back communication patterns
  // MPI_Barrier(MPI_COMM_WORLD);
  // for(int nt = 0; nt < mpi_size_; ++nt)
  //   {
  //     if(mpi_rank_ == nt)
  //       {
  // 	printf("[%d]: back\n", mpi_rank_);
  // 	printf("my_upspins.size(): %d\n", my_upspins_.size());
  // 	printf("sum_n_upspins_i_send_back_: %d\n", sum_n_upspins_i_send_back_);
  // 	printf("sum_n_upspins_i_recv_back_: %d\n", sum_n_upspins_i_recv_back_);

  // 	for(int k = 0; k < mpi_size_; ++k)
  // 	  {
  // 	    printf("-> %d (up): %d\n", k, n_upspins_i_send_back_[k]);
  // 	    int from = n_upspins_i_send_back_offsets_[k];
  // 	    int to = (k == mpi_size_-1) ? sum_n_upspins_i_send_back_ :
  // n_upspins_i_send_back_offsets_[k+1]; 	    for (int l = from; l < to;
  // ++l) 	      printf("  %s; %s\n", PrintSpinhalf(n_sites_,
  // upspins_i_send_back_[l]).c_str(), 		     PrintSpinhalf(n_sites_,
  // downspins_i_send_back_[l]).c_str());
  // 	  }

  // 	for(int k = 0; k < mpi_size_; ++k)
  // 	  {
  // 	    printf("<- %d (up): %d\n", k, n_upspins_i_recv_back_[k]);
  // 	    int from = n_upspins_i_recv_back_offsets_[k];
  // 	    int to = (k == mpi_size_-1) ? sum_n_upspins_i_recv_back_ :
  // n_upspins_i_recv_back_offsets_[k+1]; 	    for (int l = from; l < to;
  // ++l) 	      printf("  %s \n", PrintSpinhalf(n_sites_,
  // upspins_i_recv_back_[l]).c_str());
  // 	  }
  // 	printf("\n");
  //       }
  //     MPI_Barrier(MPI_COMM_WORLD);
  //   }
  // MPI_Barrier(MPI_COMM_WORLD);

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
}

template class HubbardModelMPI<double, uint32>;
template class HubbardModelMPI<complex, uint32>;
template class HubbardModelMPI<double, uint64>;
template class HubbardModelMPI<complex, uint64>;

} // namespace hydra
