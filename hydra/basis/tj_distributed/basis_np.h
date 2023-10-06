#pragma once
#ifdef HYDRA_USE_MPI

#include "extern/gsl/span"

#include <hydra/bits/bitops.h>
#include <hydra/combinatorics/lin_table.h>
#include <hydra/common.h>
#include <hydra/parallel/mpi/communicator.h>
#include <hydra/random/hash_functions.h>

namespace hydra::basis::tj_distributed {

template <typename bit_t> class BasisNp {
public:
  using bit_type = bit_t;

  BasisNp() = default;
  BasisNp(int64_t n_sites, int64_t n_up, int64_t n_dn);

  int64_t n_sites() const;
  int64_t n_up() const;
  int64_t n_dn() const;
  static constexpr bool np_conserved() { return true; }

  int64_t size() const;
  int64_t size_local() const;

  int64_t size_local_transpose() const;
  int64_t size_max() const;
  int64_t size_min() const;

  inline int rank(bit_t spins) const {
    return (int)random::hash_fnv1(spins) % mpi_size_;
  };
  mpi::Communicator transpose_communicator(bool reverse) const;

  bool operator==(BasisNp const &rhs) const;
  bool operator!=(BasisNp const &rhs) const;

private:
  int64_t n_sites_;
  int64_t n_up_;
  int64_t n_dn_;

  int64_t size_;
  int64_t size_local_;
  int64_t size_local_transpose_;
  int64_t size_max_;
  int64_t size_min_;

  int mpi_rank_;
  int mpi_size_;
  bit_t sitesmask_;

  mpi::Communicator transpose_communicator_;
  mpi::Communicator transpose_communicator_r_;

  std::vector<bit_t> my_ups_;
  std::vector<gsl::span<bit_t>> my_dns_for_ups_;
  std::vector<bit_t> my_dns_for_ups_storage_;

  std::vector<bit_t> my_dns_;
  std::vector<gsl::span<bit_t>> my_ups_for_dns_;
  std::vector<bit_t> my_ups_for_dns_storage_;

  combinatorics::LinTable<bit_t> lintable_ups_;
  combinatorics::LinTable<bit_t> lintable_dns_;
  combinatorics::LinTable<bit_t> lintable_upsc_;
  combinatorics::LinTable<bit_t> lintable_dnsc_;

public:
  std::vector<bit_t> const &my_ups() const;
  std::vector<gsl::span<bit_t>> const &my_dns_for_ups() const;
  std::vector<bit_t> const &my_dns() const;
  std::vector<gsl::span<bit_t>> const &my_ups_for_dns() const;
};

} // namespace hydra::basis::tj_distributed
#endif
