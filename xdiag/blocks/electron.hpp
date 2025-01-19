#pragma once

#include <optional>

#include <xdiag/common.hpp>

#include <xdiag/basis/electron/basis_electron.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

class ElectronIterator;

class Electron {
public:
  using basis_t = basis::BasisElectron;
  using iterator_t = ElectronIterator;

  XDIAG_API Electron() = default;
  XDIAG_API Electron(int64_t n_sites);
  XDIAG_API Electron(int64_t n_sites, int64_t nup, int64_t ndn);
  XDIAG_API Electron(int64_t n_sites, Representation const &irrep);
  XDIAG_API Electron(int64_t n_sites, int64_t nup, int64_t ndn,
                     Representation const &irrep);

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;
  XDIAG_API int64_t index(ProductState const &pstate) const;
  XDIAG_API int64_t dim() const;
  XDIAG_API int64_t size() const;

  XDIAG_API bool operator==(Electron const &rhs) const;
  XDIAG_API bool operator!=(Electron const &rhs) const;

  int64_t n_sites() const;
  std::optional<int64_t> n_up() const;
  std::optional<int64_t> n_dn() const;
  std::optional<Representation> const &irrep() const;
  bool isreal() const;
  basis_t const &basis() const;

private:
  int64_t n_sites_;
  std::optional<int64_t> n_up_;
  std::optional<int64_t> n_dn_;
  std::optional<Representation> irrep_;
  std::shared_ptr<basis_t> basis_;
  int64_t size_;
};

XDIAG_API bool isreal(Electron const &block);
XDIAG_API std::ostream &operator<<(std::ostream &out, Electron const &block);
XDIAG_API std::string to_string(Electron const &block);

class ElectronIterator {
public:
  ElectronIterator(Electron const &block, bool begin);
  ElectronIterator &operator++();
  ProductState const &operator*() const;
  bool operator!=(ElectronIterator const &rhs) const;

private:
  int64_t n_sites_;
  mutable ProductState pstate_;
  basis::BasisElectronIterator it_;
};

} // namespace xdiag
