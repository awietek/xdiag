#pragma once
#ifdef XDIAG_USE_MPI
#include <xdiag/basis/spinhalf_distributed/basis_spinhalf_distributed.hpp>
#include <xdiag/common.hpp>
#include <xdiag/states/product_state.hpp>

namespace xdiag {

class SpinhalfDistributedIterator;

class SpinhalfDistributed {
public:
  using basis_t = basis::BasisSpinhalfDistributed;
  using iterator_t = SpinhalfDistributedIterator;

  SpinhalfDistributed() = default;
  SpinhalfDistributed(int64_t n_sites, int64_t nup);

  int64_t n_sites() const;
  int64_t n_up() const;

  int64_t dim() const;
  int64_t size() const;
  int64_t size_max() const;
  int64_t size_min() const;
  iterator_t begin() const;
  iterator_t end() const;
  int64_t index(ProductState const &pstate) const;
  bool isreal(double precision = 1e-12) const;

  bool operator==(SpinhalfDistributed const &rhs) const;
  bool operator!=(SpinhalfDistributed const &rhs) const;

  basis_t const &basis() const;

private:
  int64_t n_sites_;
  int64_t n_up_;

  std::shared_ptr<basis_t> basis_;
  int64_t dim_;
  int64_t size_;
};
std::ostream &operator<<(std::ostream &out, SpinhalfDistributed const &block);
std::string to_string(SpinhalfDistributed const &block);

class SpinhalfDistributedIterator {
public:
  SpinhalfDistributedIterator(SpinhalfDistributed const &block, bool begin);
  SpinhalfDistributedIterator &operator++();
  ProductState const &operator*() const;
  bool operator!=(SpinhalfDistributedIterator const &rhs) const;

private:
  int64_t n_sites_;
  mutable ProductState pstate_;
  basis::BasisSpinhalfDistributedIterator it_;
};

} // namespace xdiag
#endif
