#pragma once

#include <hydra/common.h>

#include <hydra/combinatorics/subsets.h>

namespace hydra::indexing::electron {

template <typename bit_t> class IndexingNoNp {
public:  
  IndexingNoNp() = default;
  IndexingNoNp(int n_sites);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }

  inline idx_t size_ups() const { return size_ups_; }
  inline idx_t size_dns() const { return size_dns_; }

  inline idx_t size() const { return size_; }
  inline idx_t index_ups(bit_t ups) const { return (idx_t)ups; }
  inline idx_t index_dns(bit_t dns) const { return (idx_t)dns; }
  inline idx_t index(bit_t ups, bit_t dns) const {
    return index_ups(ups) * size_dns_ + index_dns(dns);
  }

  inline combinatorics::Subsets<bit_t> states_ups() const {
    return combinatorics::Subsets<bit_t>(n_sites_);
  }

  inline combinatorics::Subsets<bit_t> states_dns() const {
    return combinatorics::Subsets<bit_t>(n_sites_);
  }

private:
  int n_sites_;
  int n_up_;
  int n_dn_;

  idx_t size_ups_;
  idx_t size_dns_;
  idx_t size_;
};

} // namespace hydra::indexing::electron
