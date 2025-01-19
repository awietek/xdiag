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
  XDIAG_API Spinhalf(int64_t n_sites);
  XDIAG_API Spinhalf(int64_t n_sites, int64_t n_up);
  XDIAG_API Spinhalf(int64_t n_sites, Representation const &irrep);
  XDIAG_API Spinhalf(int64_t n_sites, int64_t n_up,
                     Representation const &irrep);

  // Constructors with sublattice coding
  XDIAG_API Spinhalf(int64_t n_sites, Representation const &irrep,
                     int64_t n_sublat);
  XDIAG_API Spinhalf(int64_t n_sites, int64_t n_up, Representation const &irrep,
                     int64_t n_sublat);

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;
  XDIAG_API int64_t index(ProductState const &pstate) const;
  XDIAG_API int64_t dim() const;
  XDIAG_API int64_t size() const;

  XDIAG_API bool operator==(Spinhalf const &rhs) const;
  XDIAG_API bool operator!=(Spinhalf const &rhs) const;

  int64_t n_sites() const;
  std::optional<int64_t> n_up() const;
  std::optional<Representation> const &irrep() const;
  bool isreal() const;
  basis_t const &basis() const;

private:
  int64_t n_sites_;
  std::optional<int64_t> n_up_;
  std::optional<Representation> irrep_;
  std::shared_ptr<basis_t> basis_;
  int64_t size_;
};

XDIAG_API bool isreal(Spinhalf const &block);
XDIAG_API std::ostream &operator<<(std::ostream &out, Spinhalf const &block);
XDIAG_API std::string to_string(Spinhalf const &block);

class SpinhalfIterator {
public:
  SpinhalfIterator(Spinhalf const &block, bool begin);
  SpinhalfIterator &operator++();
  ProductState const &operator*() const;
  bool operator!=(SpinhalfIterator const &rhs) const;

private:
  int64_t n_sites_;
  mutable ProductState pstate_;
  basis::BasisSpinhalfIterator it_;
};
} // namespace xdiag
