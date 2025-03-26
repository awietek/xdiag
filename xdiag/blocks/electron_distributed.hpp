#pragma once
#ifdef XDIAG_USE_MPI
#include <optional>

#include <xdiag/basis/electron_distributed/basis_electron_distributed.hpp>
#include <xdiag/common.hpp>
#include <xdiag/states/product_state.hpp>

namespace xdiag {

class ElectronDistributedIterator;

class ElectronDistributed {
public:
  using basis_t = basis::BasisElectronDistributed;
  using iterator_t = ElectronDistributedIterator;

  XDIAG_API ElectronDistributed() = default;
  XDIAG_API ElectronDistributed(int64_t nsites, int64_t nup, int64_t ndn,
                                std::string backend = "auto");

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;
  XDIAG_API int64_t index(ProductState const &pstate) const;
  XDIAG_API int64_t dim() const;
  XDIAG_API int64_t size() const;
  XDIAG_API int64_t size_max() const;
  XDIAG_API int64_t size_min() const;

  XDIAG_API bool operator==(ElectronDistributed const &rhs) const;
  XDIAG_API bool operator!=(ElectronDistributed const &rhs) const;

  XDIAG_API int64_t nsites() const;
  XDIAG_API bool isreal() const;

  std::string backend() const;
  std::optional<int64_t> nup() const;
  std::optional<int64_t> ndn() const;
  basis_t const &basis() const;

private:
  int64_t nsites_;
  std::string backend_;
  std::optional<int64_t> nup_;
  std::optional<int64_t> ndn_;
  std::shared_ptr<basis_t> basis_;
  int64_t dim_;
  int64_t size_;
};

XDIAG_API int64_t index(ElectronDistributed const &block,
                        ProductState const &pstate);
XDIAG_API int64_t nsites(ElectronDistributed const &block);
XDIAG_API int64_t dim(ElectronDistributed const &block);
XDIAG_API int64_t size(ElectronDistributed const &block);
XDIAG_API bool isreal(ElectronDistributed const &block);
XDIAG_API std::ostream &operator<<(std::ostream &out,
                                   ElectronDistributed const &block);
XDIAG_API std::string to_string(ElectronDistributed const &block);

class ElectronDistributedIterator {
public:
  ElectronDistributedIterator(ElectronDistributed const &block, bool begin);
  XDIAG_API ElectronDistributedIterator &operator++();
  XDIAG_API ProductState const &operator*() const;
  XDIAG_API bool operator!=(ElectronDistributedIterator const &rhs) const;

private:
  int64_t nsites_;
  mutable ProductState pstate_;
  basis::BasisElectronDistributedIterator it_;
};

} // namespace xdiag

#endif
