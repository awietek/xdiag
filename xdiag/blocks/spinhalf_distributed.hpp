#pragma once
#ifdef XDIAG_USE_MPI
#include <optional>

#include <xdiag/basis/spinhalf_distributed/basis_spinhalf_distributed.hpp>
#include <xdiag/common.hpp>
#include <xdiag/states/product_state.hpp>

namespace xdiag {

class SpinhalfDistributedIterator;

class SpinhalfDistributed {
public:
  using basis_t = basis::BasisSpinhalfDistributed;
  using iterator_t = SpinhalfDistributedIterator;

  XDIAG_API SpinhalfDistributed() = default;
  XDIAG_API SpinhalfDistributed(int64_t nsites, int64_t nup,
                                std::string backend = "auto");

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;
  XDIAG_API int64_t index(ProductState const &pstate) const;
  XDIAG_API int64_t dim() const;
  XDIAG_API int64_t size() const;
  XDIAG_API int64_t size_max() const;
  XDIAG_API int64_t size_min() const;

  XDIAG_API bool operator==(SpinhalfDistributed const &rhs) const;
  XDIAG_API bool operator!=(SpinhalfDistributed const &rhs) const;

  XDIAG_API int64_t nsites() const;
  XDIAG_API bool isreal(double precision = 1e-12) const;

  std::string backend() const;
  std::optional<int64_t> nup() const;
  basis_t const &basis() const;

private:
  int64_t nsites_;
  std::string backend_;
  std::optional<int64_t> nup_;
  std::shared_ptr<basis_t> basis_;
  int64_t dim_;
  int64_t size_;
};

XDIAG_API int64_t index(SpinhalfDistributed const &block,
                        ProductState const &pstate);
XDIAG_API int64_t nsites(SpinhalfDistributed const &block);
XDIAG_API int64_t dim(SpinhalfDistributed const &block);
XDIAG_API int64_t size(SpinhalfDistributed const &block);
XDIAG_API bool isreal(SpinhalfDistributed const &block);
XDIAG_API std::ostream &operator<<(std::ostream &out,
                                   SpinhalfDistributed const &block);
XDIAG_API std::string to_string(SpinhalfDistributed const &block);

class SpinhalfDistributedIterator {
public:
  SpinhalfDistributedIterator(SpinhalfDistributed const &block, bool begin);
  XDIAG_API SpinhalfDistributedIterator &operator++();
  XDIAG_API ProductState const &operator*() const;
  XDIAG_API bool operator!=(SpinhalfDistributedIterator const &rhs) const;

private:
  int64_t nsites_;
  mutable ProductState pstate_;
  basis::BasisSpinhalfDistributedIterator it_;
};

} // namespace xdiag
#endif
