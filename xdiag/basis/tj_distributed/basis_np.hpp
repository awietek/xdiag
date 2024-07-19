#pragma once
#ifdef XDIAG_USE_MPI

#include <unordered_map>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/common.hpp>
#include <xdiag/extern/gsl/span>
#include <xdiag/parallel/mpi/communicator.hpp>
#include <xdiag/random/hash_functions.hpp>

namespace xdiag::basis::tj_distributed {

template <typename bit_tt> class BasisNp {
public:
  using bit_t = bit_tt;

  BasisNp() = default;
  BasisNp(int64_t n_sites, int64_t n_up, int64_t n_dn);

  int64_t n_sites() const;
  int64_t n_up() const;
  int64_t n_dn() const;
  static constexpr bool np_conserved() { return true; }

  int64_t index(bit_t up, bit_t dn) const;

  int64_t dim() const;
  int64_t size() const;
  int64_t size_transpose() const;
  int64_t size_max() const;
  int64_t size_min() const;

  bool operator==(BasisNp const &rhs) const;
  bool operator!=(BasisNp const &rhs) const;

private:
  int64_t n_sites_;
  int64_t n_up_;
  int64_t n_dn_;

  combinatorics::LinTable<bit_t> lintable_dncs_;
  combinatorics::LinTable<bit_t> lintable_upcs_;

  int64_t dim_;
  int64_t size_;
  int64_t size_transpose_;
  int64_t size_max_;
  int64_t size_min_;

  int mpi_rank_;
  int mpi_size_;
  bit_t sitesmask_;

  mpi::Communicator transpose_communicator_;
  mpi::Communicator transpose_communicator_r_;
  std::vector<int64_t> transpose_permutation_;
  std::vector<int64_t> transpose_permutation_r_;

  std::vector<bit_t> my_ups_;
  std::unordered_map<bit_t, int64_t> my_ups_offset_;
  std::vector<gsl::span<bit_t>> my_dns_for_ups_;
  std::vector<bit_t> my_dns_for_ups_storage_;

  std::vector<bit_t> my_dns_;
  std::unordered_map<bit_t, int64_t> my_dns_offset_;
  std::vector<gsl::span<bit_t>> my_ups_for_dns_;
  std::vector<bit_t> my_ups_for_dns_storage_;

public:
  std::vector<bit_t> const &my_ups() const;
  int64_t my_ups_offset(bit_t ups) const;
  gsl::span<bit_t> my_dns_for_ups(int64_t idx_ups) const;
  inline bit_t my_dns_for_ups_storage(int64_t idx) const {
    return my_dns_for_ups_storage_[idx];
  }

  std::vector<bit_t> const &my_dns() const;
  int64_t my_dns_offset(bit_t dns) const;
  gsl::span<bit_t> my_ups_for_dns(int64_t idx_dns) const;
  inline bit_t my_ups_for_dns_storage(int64_t idx) const {
    return my_ups_for_dns_storage_[idx];
  }

  inline int rank(bit_t spins) const { // mpi ranks are ints
    return (int)(random::hash_div3(spins) % mpi_size_);
  };
  inline int64_t index_dncs(bit_t dncs) const {
    return lintable_dncs_.index(dncs);
  }
  inline int64_t index_upcs(bit_t upcs) const {
    return lintable_upcs_.index(upcs);
  }

  // transforms a vector in up/dn order to dn/up order
  // if no "out_vec" is given result of transpose is stored
  // in send_buffer of mpi::buffer
  // recv_buffer is filled with zeros
  template <typename coeff_t>
  void transpose(const coeff_t *in_vec, coeff_t *out_vec = nullptr) const;

  // transforms a vector in dn/up order to up/dn order
  // if no "out_vec" is given result of transpose is stored
  // in send_buffer of mpi::buffer
  // recv_buffer is filled with zeros
  template <typename coeff_t>
  void transpose_r(coeff_t const *in_vec, coeff_t *out_vec = nullptr) const;
};

} // namespace xdiag::basis::tj_distributed
#endif
