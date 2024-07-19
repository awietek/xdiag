#pragma once
#ifdef XDIAG_USE_MPI

#include <xdiag/extern/gsl/span>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/common.hpp>
#include <xdiag/parallel/mpi/comm_pattern.hpp>
#include <xdiag/parallel/mpi/communicator.hpp>
#include <xdiag/random/hash_functions.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <typename bit_tt> class BasisSz {
public:
  using bit_t = bit_tt;

  BasisSz() = default;
  BasisSz(int64_t n_sites, int64_t n_up);

  int64_t n_sites() const;
  int64_t n_up() const;

  int64_t n_prefix_bits() const;
  int64_t n_postfix_bits() const;

  int64_t dim() const;
  int64_t size() const;
  int64_t size_transpose() const;
  int64_t size_max() const;
  int64_t size_min() const;

  std::vector<bit_t> const &prefixes() const;
  int64_t prefix_begin(bit_t prefix) const;
  combinatorics::LinTable<bit_t> const &postfix_lintable(bit_t prefix) const;
  std::vector<bit_t> const &postfix_states(bit_t prefix) const;

  std::vector<bit_t> const &postfixes() const;
  int64_t postfix_begin(bit_t postfix) const;
  combinatorics::LinTable<bit_t> const &prefix_lintable(bit_t postfix) const;
  std::vector<bit_t> const &prefix_states(bit_t postfix) const;

  inline int rank(bit_t spins) const {
    return (int)(random::hash_div3(spins) % mpi_size_);
  };

  mpi::CommPattern &comm_pattern();
  mpi::Communicator transpose_communicator(bool reverse) const;

  bool operator==(BasisSz const &rhs) const;
  bool operator!=(BasisSz const &rhs) const;

private:
  int64_t n_sites_;
  int64_t n_up_;

  int n_prefix_bits_;
  int n_postfix_bits_;

  int64_t dim_;
  int64_t size_;
  int64_t size_transpose_;
  int64_t size_max_;
  int64_t size_min_;

  int mpi_rank_;
  int mpi_size_;

  std::vector<bit_t> prefixes_;
  std::unordered_map<bit_t, int64_t> prefix_begin_;
  std::vector<combinatorics::LinTable<bit_t>> postfix_lintables_;
  std::vector<std::vector<bit_t>> postfix_states_;

  std::vector<bit_t> postfixes_;
  std::unordered_map<bit_t, int64_t> postfix_begin_;
  std::vector<combinatorics::LinTable<bit_t>> prefix_lintables_;
  std::vector<std::vector<bit_t>> prefix_states_;

  mpi::CommPattern comm_pattern_;
  mpi::Communicator transpose_communicator_;
  mpi::Communicator transpose_communicator_reverse_;
};

} // namespace xdiag::basis::spinhalf_distributed
#endif
