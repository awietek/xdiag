#pragma once
#ifdef XDIAG_USE_MPI

#include <xdiag/extern/gsl/span>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/common.hpp>
#include <xdiag/parallel/mpi/communicator.hpp>
#include <xdiag/random/hash_functions.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <typename bit_t> class BasisSz {
public:
  using bit_type = bit_t;

  BasisSz() = default;
  BasisSz(int64_t n_sites, int64_t n_up, int64_t n_dn);

  int64_t n_sites() const;
  int64_t n_up() const;
  static constexpr bool sz_conserved() { return true; }

  int64_t size() const;
  int64_t size_local() const;
  int64_t size_local_transpose() const;
  int64_t size_max() const;
  int64_t size_min() const;

  inline int rank(bit_t spins) const {
    return (int)random::hash_fnv1(spins) % mpi_size_;
  };
  mpi::Communicator transpose_communicator(bool reverse) const;

  bool operator==(BasisSz const &rhs) const;
  bool operator!=(BasisSz const &rhs) const;

private:
  int64_t n_sites_;
  int64_t n_up_;

  int64_t size_;
  int64_t size_local_;
  int64_t size_local_transpose_;
  int64_t size_max_;
  int64_t size_min_;

  int mpi_rank_;
  int mpi_size_;

  int n_prefix_bits_;
  int n_postfix_bits_;

  std::vector<bit_t> prefixes_;
  std::unordered_map<bit_t, int64_t> prefix_begin_;
  std::vector<indexing::LinTable<bit_t>> postfix_lintables_;
  std::vector<std::vector<bit_t>> postfix_states_;

  std::vector<bit_t> postfixes_;
  std::unordered_map<bit_t, int64_t> postfix_begin_;
  std::vector<indexing::LinTable<bit_t>> prefix_lintables_;
  std::vector<std::vector<bit_t>> prefix_states_;

public:

  
};

} // namespace xdiag::basis::tj_distributed
#endif
