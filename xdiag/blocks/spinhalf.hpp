#pragma once

#include <optional>

#include <xdiag/common.hpp>

#include <xdiag/basis/spinhalf/basis_spinhalf.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

class SpinhalfIterator;

class Spinhalf {
public:
  using basis_t = basis::BasisSpinhalf;
  using iterator_t = SpinhalfIterator;

  XDIAG_API Spinhalf() = default;
  XDIAG_API Spinhalf(int64_t nsites, std::string backend = "auto");
  XDIAG_API Spinhalf(int64_t nsites, int64_t nup, std::string backend = "auto");
  XDIAG_API Spinhalf(int64_t nsites, Representation const &irrep,
                     std::string backend = "auto");
  XDIAG_API Spinhalf(int64_t nsites, int64_t nup, Representation const &irrep,
                     std::string backend = "auto");

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;
  XDIAG_API int64_t index(ProductState const &pstate) const;
  XDIAG_API int64_t dim() const;
  XDIAG_API int64_t size() const;

  XDIAG_API bool operator==(Spinhalf const &rhs) const;
  XDIAG_API bool operator!=(Spinhalf const &rhs) const;

  XDIAG_API int64_t nsites() const;
  XDIAG_API bool isreal() const;

  std::string backend() const;
  std::optional<int64_t> nup() const;
  std::optional<Representation> const &irrep() const;
  basis_t const &basis() const;
private:
  int64_t nsites_;
  std::string backend_;
  std::optional<int64_t> nup_;
  std::optional<Representation> irrep_;
  std::shared_ptr<basis_t> basis_;
  int64_t size_;
};

XDIAG_API int64_t index(Spinhalf const &block, ProductState const &pstate);
XDIAG_API int64_t nsites(Spinhalf const &block);
XDIAG_API int64_t dim(Spinhalf const &block);
XDIAG_API int64_t size(Spinhalf const &block);
XDIAG_API bool isreal(Spinhalf const &block);
XDIAG_API std::ostream &operator<<(std::ostream &out, Spinhalf const &block);
XDIAG_API std::string to_string(Spinhalf const &block);

class SpinhalfIterator {
public:
  SpinhalfIterator(Spinhalf const &block, bool begin);
  XDIAG_API SpinhalfIterator &operator++();
  XDIAG_API ProductState const &operator*() const;
  XDIAG_API bool operator!=(SpinhalfIterator const &rhs) const;

private:
  int64_t nsites_;
  mutable ProductState pstate_;
  basis::BasisSpinhalfIterator it_;
};
} // namespace xdiag
