#pragma once

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>

namespace hydra::indexing {

template <typename bit_t> class ElectronIndexing {
public:
  ElectronIndexing() = default;
  ElectronIndexing(int n_sites, int n_up, int n_dn);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }

  inline idx_t size_ups() const { return size_ups_; }
  inline idx_t size_dns() const { return size_dns_; }

  inline idx_t size() const { return size_; }
  inline idx_t index_ups(bit_t ups) const {
    return lintable_ups_.index(ups);
  }
  inline idx_t index_dns(bit_t dns) const {
    return lintable_dns_.index(dns);
  }

private:
  int n_sites_;
  int n_up_;
  int n_dn_;

  idx_t size_ups_;
  idx_t size_dns_;
  idx_t size_;

  indexing::LinTable<bit_t> lintable_ups_;
  indexing::LinTable<bit_t> lintable_dns_;
};

} // namespace hydra::indexing
