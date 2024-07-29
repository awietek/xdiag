#pragma once
#ifdef XDIAG_USE_MPI
#include <xdiag/basis/tj_distributed/basis_tj_distributed.hpp>
#include <xdiag/common.hpp>
#include <xdiag/states/product_state.hpp>

namespace xdiag {

class tJDistributedIterator;

class tJDistributed {
public:
  using basis_t = basis::BasistJDistributed;
  using iterator_t = tJDistributedIterator;
  tJDistributed() = default;
  tJDistributed(int64_t n_sites, int64_t n_up, int64_t n_dn);

  int64_t n_sites() const;
  int64_t n_up() const;
  int64_t n_dn() const;

  int64_t dim() const;
  int64_t size() const;
  int64_t size_max() const;
  int64_t size_min() const;
  iterator_t begin() const;
  iterator_t end() const;
  int64_t index(ProductState const &pstate) const;
  bool isreal(double precision = 1e-12) const;

  bool operator==(tJDistributed const &rhs) const;
  bool operator!=(tJDistributed const &rhs) const;
  basis_t const &basis() const;

private:
  int64_t n_sites_;
  int64_t n_up_;
  int64_t n_dn_;

  std::shared_ptr<basis_t> basis_;
  int64_t dim_;
  int64_t size_;
};

std::ostream &operator<<(std::ostream &out, tJDistributed const &block);
std::string to_string(tJDistributed const &block);

class tJDistributedIterator {
public:
  tJDistributedIterator(tJDistributed const &block, bool begin);
  tJDistributedIterator &operator++();
  ProductState const &operator*() const;
  bool operator!=(tJDistributedIterator const &rhs) const;

private:
  int64_t n_sites_;
  mutable ProductState pstate_;
  basis::BasistJDistributedIterator it_;
};

} // namespace xdiag

#endif
