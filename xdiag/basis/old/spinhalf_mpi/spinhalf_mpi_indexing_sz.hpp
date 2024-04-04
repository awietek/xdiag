#pragma once

#include <mpi.h>
#include <unordered_map>

#include <xdiag/common.hpp>
#include <lila/external/gsl/span>

#include <xdiag/combinatorics/combinations.hpp>

#include <xdiag/indexing/lintable.hpp>

namespace xdiag::indexing {
template <typename bit_t> class SpinhalfMPIIndexingSz {
public:
  SpinhalfMPIIndexingSz() = default;
  SpinhalfMPIIndexingSz(int n_sites, int nup);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_prefix_bits() const { return n_prefix_bits_; }
  inline int n_postfix_bits() const { return n_postfix_bits_; }
  inline int64_t size() const { return size_; }
  inline int64_t size_transpose() const { return size_transpose_; }
  inline int64_t size_max() const { return size_max_; }
  inline int64_t dim() const { return dim_; }

  int process(bit_t prepostfix) const;
  inline int mpi_rank() const { return mpi_rank_; }
  inline int mpi_size() const { return mpi_size_; }

  gsl::span<bit_t const> prefixes() const;
  int64_t prefix_begin(bit_t prefix) const {
    // for (auto [p, b] : prefix_begin_)
    // 	lila::Log("pi: {} p: {} b: {}", prefix, p, b);
    if (prefix_begin_.find(prefix) != prefix_begin_.end()) {
      return prefix_begin_.at(prefix);
    } else {
      return invalid_index;
    }
  }
  gsl::span<bit_t const> postfixes(bit_t prefix) const;
  LinTable<bit_t> const &postfix_indexing(bit_t prefix) const;

  gsl::span<bit_t const> postfixes() const;
  int64_t postfix_begin(bit_t postfix) const {
    // for (auto [p, b] : postfix_begin_)
    //   lila::Log("#{} pi: {} p: {} b: {}, presize: {}", mpi_rank_, postfix, p,
    //   b, prefixes(postfix).size());
    if (postfix_begin_.find(postfix) != postfix_begin_.end()) {
      return postfix_begin_.at(postfix);
    } else {
      return invalid_index;
    }
  }
  gsl::span<bit_t const> prefixes(bit_t postfix) const;
  LinTable<bit_t> const &prefix_indexing(bit_t postfix) const;

  LinTable<bit_t> const &empty_lintable() const { return empty_lintable_; }

private:
  int n_sites_;
  int n_up_;

  int n_prefix_bits_;
  int n_postfix_bits_;

  indexing::LinTable<bit_t> empty_lintable_;

  std::vector<bit_t> prefixes_;
  std::unordered_map<bit_t, int64_t> prefix_begin_;
  std::vector<indexing::LinTable<bit_t>> postfix_lintables_;
  std::vector<std::vector<bit_t>> postfix_states_;

  std::vector<bit_t> postfixes_;
  std::unordered_map<bit_t, int64_t> postfix_begin_;
  std::vector<indexing::LinTable<bit_t>> prefix_lintables_;
  std::vector<std::vector<bit_t>> prefix_states_;

  int64_t size_;
  int64_t size_transpose_;
  int64_t size_max_;
  int64_t dim_;

  int mpi_rank_;
  int mpi_size_;
};

} // namespace xdiag::indexing
