#include "reduced_density_matrix.h"

#include <hydra/indexing/index_table.h>
#include <hydra/mpi/communicator.h>

namespace hydra {

template <class coeff_t>
lila::Matrix<coeff_t> ReducedDensityMatrix(TJModel<coeff_t> const &model,
                                           int n_sites_A, qn_tj qn_A,
                                           lila::Vector<coeff_t> const &wf) {
  using model_t = TJModel<coeff_t>;
  using idx_t = typename model_t::idx_t;
  using basis_t = typename model_t::basis_t;

  auto qn_tot = model.qn();
  auto qn_B = qn_tot - qn_A;
  int n_sites_tot = model.n_sites();
  int n_sites_B = n_sites_tot - n_sites_A;

  assert(n_sites_tot >= 0);
  assert(n_sites_A >= 0);
  assert(n_sites_B >= 0);

  assert(valid(qn_tot, n_sites_tot));
  assert(valid(qn_A, n_sites_A));
  assert(valid(qn_B, n_sites_B));

  auto basis_A = basis_t(n_sites_A, qn_A);
  auto basis_B = basis_t(n_sites_B, qn_B);
  auto basis_tot = basis_t(n_sites_tot, qn_tot);
  IndexTable<basis_t, idx_t> indexing(basis_tot);

  // Compute the partial trace over B
  auto dim = basis_A.size();
  auto rho_A = lila::Zeros<coeff_t>(dim, dim);

  idx_t i_A1 = 0;
  for (auto state_A1 : basis_A) {
    idx_t i_A2 = 0;
    for (auto state_A2 : basis_A) {
      for (auto state_B : basis_B) {
        auto state_B_shifted = (state_B << n_sites_A);
        auto state_1 = state_A1 | state_B_shifted;
        auto state_2 = state_A2 | state_B_shifted;
        auto idx_1 = indexing.index(state_1);
        auto idx_2 = indexing.index(state_2);

        rho_A(i_A1, i_A2) += lila::conj(wf(idx_1)) * wf(idx_2);
      }
      ++i_A2;
    }
    ++i_A1;
  }

  return rho_A;
}

template <class coeff_t>
lila::Matrix<coeff_t> ReducedDensityMatrix(TJModelMPI<coeff_t> const &model,
                                           int n_sites_A, qn_tj qn_A,
                                           lila::VectorMPI<coeff_t> const &wf) {

  using model_t = TJModel<coeff_t>;
  using bit_t = typename model_t::idx_t;
  using idx_t = typename model_t::idx_t;
  using basis_t = typename model_t::basis_t;
  using state_t = typename model_t::state_t;

  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Get bases for A/B subsystem
  auto qn_tot = model.qn();
  auto qn_B = qn_tot - qn_A;
  int n_sites_tot = model.n_sites();
  int n_sites_B = n_sites_tot - n_sites_A;

  assert(n_sites_tot >= 0);
  assert(n_sites_A >= 0);
  assert(n_sites_B >= 0);

  assert(valid(qn_tot, n_sites_tot));
  assert(valid(qn_A, n_sites_A));
  assert(valid(qn_B, n_sites_B));

  auto basis_A = basis_t(n_sites_A, qn_A);
  auto basis_B = basis_t(n_sites_B, qn_B);

  // Determine communication pattern
  auto process_no = [&mpi_size, &n_sites_A](state_t const &state) {
    auto state_B = state >> n_sites_A;
    return hash(state_B) % mpi_size;
  };
  std::vector<idx_t> n_states_i_send(mpi_size, 0);
  model.transverse(
      [process_no, &n_states_i_send](state_t const &state, idx_t &idx) {
        auto proc = process_no(state);
        ++n_states_i_send[proc];
      });
  mpi::Communicator<idx_t> comm_ups(n_states_i_send);
  mpi::Communicator<idx_t> comm_dns(n_states_i_send);
  mpi::Communicator<idx_t> comm_coeff(n_states_i_send);

  // Fill send buffers
  std::vector<bit_t> ups_buffer_send(comm_ups.send_buffer_size());
  std::vector<bit_t> dns_buffer_send(comm_dns.send_buffer_size());
  std::vector<coeff_t> coeff_buffer_send(comm_coeff.send_buffer_size());

  model.transverse([process_no, &ups_buffer_send, &dns_buffer_send,
                    &coeff_buffer_send, &comm_ups, &comm_dns, &comm_coeff,
                    &wf](state_t const &state, idx_t &idx) {
    auto proc = process_no(state);
    comm_ups.add_to_send_buffer(proc, state.ups, ups_buffer_send.data());
    comm_dns.add_to_send_buffer(proc, state.dns, dns_buffer_send.data());
    comm_coeff.add_to_send_buffer(proc, wf(idx), coeff_buffer_send.data());
  });

  // Send states around
  std::vector<bit_t> ups_buffer_recv(comm_ups.recv_buffer_size());
  std::vector<bit_t> dns_buffer_recv(comm_dns.recv_buffer_size());
  std::vector<coeff_t> coeff_buffer_recv(comm_coeff.recv_buffer_size());
  comm_ups.all_to_all(ups_buffer_send.data(), ups_buffer_recv.data());
  comm_dns.all_to_all(dns_buffer_send.data(), dns_buffer_recv.data());
  comm_coeff.all_to_all(coeff_buffer_send.data(), coeff_buffer_recv.data());
  ups_buffer_send.resize(0);
  dns_buffer_send.resize(0);
  coeff_buffer_send.resize(0);
  ups_buffer_send.shrink_to_fit();
  dns_buffer_send.shrink_to_fit();
  coeff_buffer_send.shrink_to_fit();

  // create vector of states and coeffs
  std::vector<coeff_t> &coeffs = coeff_buffer_recv;
  std::vector<state_t> states;
  for (unsigned i = 0; i < ups_buffer_recv.size(); ++i) {
    state_t state{ups_buffer_recv[i], dns_buffer_recv[i]}; 
    states.push_back(state);
    assert(process_no(state) == mpi_rank);
  }
  ups_buffer_recv.resize(0);
  dns_buffer_recv.resize(0);
  ups_buffer_recv.shrink_to_fit();
  dns_buffer_recv.shrink_to_fit();

  // Sort states and coeffs
  std::vector<size_t> idx(states.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&states](size_t i1, size_t i2) {
    return states[i1] < states[i2];
  });

  std::vector<state_t> states_sorted;
  std::vector<coeff_t> coeffs_sorted;
  states_sorted.reserve(states.size());
  coeffs_sorted.reserve(states.size());
  for (auto j : idx) {
    states_sorted.emplace_back(states[j]);
    coeffs_sorted.emplace_back(coeffs[j]);
  }
  states.resize(0);
  states.shrink_to_fit();
  coeffs.resize(0);
  coeffs.shrink_to_fit();

  auto index = [&states_sorted](state_t const &state) {
    auto it =
        std::lower_bound(states_sorted.begin(), states_sorted.end(), state);
    assert(std::find(states_sorted.begin(), states_sorted.end(), state) !=
           states_sorted.end());
    return std::distance(states_sorted.begin(), it);
  };

  // Compute the "partial" partial trace over B
  auto dim = basis_A.size();
  auto rho_A = lila::Zeros<coeff_t>(dim, dim);
  idx_t i_A1 = 0;
  for (auto state_A1 : basis_A) {
    idx_t i_A2 = 0;
    for (auto state_A2 : basis_A) {
      for (auto state_B : basis_B) {
        auto state_B_shifted = (state_B << n_sites_A);
        if (process_no(state_B_shifted) == mpi_rank) {
          auto state_1 = state_A1 | state_B_shifted;
          auto state_2 = state_A2 | state_B_shifted;
	  assert(process_no(state_B_shifted) == mpi_rank);
	  assert(process_no(state_1) == mpi_rank);
	  assert(process_no(state_2) == mpi_rank);
          auto idx_1 = index(state_1);
          auto idx_2 = index(state_2);
          rho_A(i_A1, i_A2) +=
              lila::conj(coeffs_sorted[idx_1]) * coeffs_sorted[idx_2];
        }
      }
      ++i_A2;
    }
    ++i_A1;
  }

  // Sum up "partial" partial traces over B
  auto rho_A2 = lila::Zeros<coeff_t>(dim, dim);
  lila::MPI_Allreduce<coeff_t>(rho_A.data(), rho_A2.data(), rho_A.size(),
                               MPI_SUM, MPI_COMM_WORLD);

  // LilaPrint(rho_A);
  return rho_A2;
}

template lila::Matrix<double>
ReducedDensityMatrix(TJModel<double> const &model, int n_sites_A, qn_tj qn_A,
                     lila::Vector<double> const &wf);

template lila::Matrix<complex>
ReducedDensityMatrix(TJModel<complex> const &model, int n_sites_A, qn_tj qn_A,
                     lila::Vector<complex> const &wf);

template lila::Matrix<double>
ReducedDensityMatrix(TJModelMPI<double> const &model, int n_sites_A, qn_tj qn_A,
                     lila::VectorMPI<double> const &wf);

template lila::Matrix<complex>
ReducedDensityMatrix(TJModelMPI<complex> const &model, int n_sites_A,
                     qn_tj qn_A, lila::VectorMPI<complex> const &wf);

} // namespace hydra
