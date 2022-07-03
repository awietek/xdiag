#pragma once

#include <unordered_map>

#include <hydra/common.h>
#include <lila/external/gsl/span>

#include <hydra/combinatorics/combinations.h>

#include <hydra/indexing/lintable.h>

namespace hydra::indexing {
template <typename bit_t> class SpinhalfMPIIndexingSz {
public:
  SpinhalfMPIIndexingSz() = default;
  SpinhalfMPIIndexingSz(int n_sites, int nup);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_prefix_bits() const { return n_prefix_bits_; }
  inline int n_postfix_bits() const { return n_postfix_bits_; }
  inline idx_t size() const { return size_; }
  inline idx_t size_transpose() const { return size_transpose_; }
  inline idx_t dim() const { return dim_; }

  int process(bit_t prepostfix) const;
  inline int mpi_rank() const { return mpi_rank_; }
  inline int mpi_size() const { return mpi_size_; }
  
  gsl::span<bit_t const> prefixes() const;
  idx_t prefix_begin(bit_t prefix) const { return prefix_begin_.at(prefix); }
  gsl::span<bit_t const> postfixes(bit_t prefix) const;
  LinTable<bit_t> const &postfix_indexing(bit_t prefix) const;

  gsl::span<bit_t const> postfixes() const;
  idx_t postfix_begin(bit_t postfix) const {
    return postfix_begin_.at(postfix);
  }
  gsl::span<bit_t const> prefixes(bit_t postfix) const;
  LinTable<bit_t> const &prefix_indexing(bit_t postfix) const;

private:
  int n_sites_;
  int n_up_;

  int n_prefix_bits_;
  int n_postfix_bits_;

  std::vector<bit_t> prefixes_;
  std::unordered_map<bit_t, idx_t> prefix_begin_;
  std::vector<indexing::LinTable<bit_t>> postfix_lintables_;
  std::vector<std::vector<bit_t>> postfix_states_;

  std::vector<bit_t> postfixes_;
  std::unordered_map<bit_t, idx_t> postfix_begin_;
  std::vector<indexing::LinTable<bit_t>> prefix_lintables_;
  std::vector<std::vector<bit_t>> prefix_states_;

  idx_t size_;
  idx_t size_transpose_;
  idx_t dim_;

  int mpi_rank_;
  int mpi_size_;
};

} // namespace hydra::indexing
