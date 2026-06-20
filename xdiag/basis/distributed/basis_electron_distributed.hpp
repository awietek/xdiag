// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <unordered_map>

#include <extern/gsl/span>

#include <xdiag/basis/basis.hpp>
#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/extract_deposit.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/mpi/communicator.hpp>
#include <xdiag/random/hash_functions.hpp>
#include <xdiag/utils/type_name.hpp>

namespace xdiag::basis {

template <typename bit_tt> class BasisElectronDistributedIterator;

template <typename bit_tt>
class BasisElectronDistributed
    : public BasisType<BasisElectronDistributed<bit_tt>> {
public:
  using bit_t = bit_tt;
  using iterator_t = BasisElectronDistributedIterator<bit_t>;
  static constexpr std::string_view type_name =
      utils::get_type_name<BasisElectronDistributed<bit_t>>();

  BasisElectronDistributed() = default;
  BasisElectronDistributed(int64_t nsites, int64_t nup, int64_t ndn);

  int64_t nsites() const;
  int64_t nup() const;
  int64_t ndn() const;
  // static constexpr bool np_conserved() { return true; }

  int64_t dim() const;
  int64_t size() const;
  int64_t size_transpose() const;
  int64_t size_max() const;
  int64_t size_min() const;
  iterator_t begin() const;
  iterator_t end() const;
  int64_t index(bit_t up, bit_t dn) const;

  bool operator==(BasisElectronDistributed const &rhs) const;
  bool operator!=(BasisElectronDistributed const &rhs) const;

private:
  int64_t nsites_;
  int64_t nup_;
  int64_t ndn_;

  combinatorics::LinTable<bit_t> lintable_dns_;
  combinatorics::LinTable<bit_t> lintable_ups_;

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
  std::vector<bit_t> my_dns_;
  std::unordered_map<bit_t, int64_t> my_dns_offset_;

  std::vector<bit_t> all_ups_;
  std::vector<bit_t> all_dns_;

public:
  std::vector<bit_t> const &my_ups() const;
  int64_t my_ups_offset(bit_t ups) const;
  std::vector<bit_t> const &all_dns() const;

  std::vector<bit_t> const &my_dns() const;
  int64_t my_dns_offset(bit_t dns) const;
  std::vector<bit_t> const &all_ups() const;

  inline int rank(bit_t spins) const { // mpi ranks are ints
    return (int)(random::hash_div3(spins) % mpi_size_);
  };
  inline int64_t index_dns(bit_t dns) const { return lintable_dns_.index(dns); }
  inline int64_t index_ups(bit_t ups) const { return lintable_ups_.index(ups); }

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

template <typename bit_tt> class BasisElectronDistributedIterator {
public:
  using bit_t = bit_tt;
  BasisElectronDistributedIterator() = default;
  BasisElectronDistributedIterator(BasisElectronDistributed<bit_t> const &basis,
                                   bool begin);
  BasisElectronDistributedIterator<bit_t> &operator++();
  std::pair<bit_t, bit_t> operator*() const;
  bool operator!=(BasisElectronDistributedIterator<bit_t> const &rhs) const;

private:
  BasisElectronDistributed<bit_t> const &basis_;
  int64_t up_idx_;
  int64_t dn_idx_;
};

} // namespace xdiag::basis
