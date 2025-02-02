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
  XDIAG_API Electron(int64_t nsites, std::string backend = "auto");
  XDIAG_API Electron(int64_t nsites, int64_t nup, int64_t ndn,
                     std::string backend = "auto");
  XDIAG_API Electron(int64_t nsites, Representation const &irrep,
                     std::string backend = "auto");
  XDIAG_API Electron(int64_t nsites, int64_t nup, int64_t ndn,
                     Representation const &irrep, std::string backend = "auto");

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;
  XDIAG_API int64_t index(ProductState const &pstate) const;
  XDIAG_API int64_t dim() const;
  XDIAG_API int64_t size() const;

  XDIAG_API bool operator==(Electron const &rhs) const;
  XDIAG_API bool operator!=(Electron const &rhs) const;

  XDIAG_API int64_t nsites() const;
  XDIAG_API bool isreal() const;

  std::string backend() const;
  std::optional<int64_t> nup() const;
  std::optional<int64_t> ndn() const;
  std::optional<Representation> const &irrep() const;
  basis_t const &basis() const;

private:
  int64_t nsites_;
  std::string backend_;
  std::optional<int64_t> nup_;
  std::optional<int64_t> ndn_;
  std::optional<Representation> irrep_;
  std::shared_ptr<basis_t> basis_;
  int64_t size_;
};

XDIAG_API int64_t index(Electron const &block, ProductState const &pstate);
XDIAG_API int64_t nsites(Electron const &block);
XDIAG_API int64_t dim(Electron const &block);
XDIAG_API int64_t size(Electron const &block);
XDIAG_API bool isreal(Electron const &block);
XDIAG_API std::ostream &operator<<(std::ostream &out, Electron const &block);
XDIAG_API std::string to_string(Electron const &block);

class ElectronIterator {
public:
  ElectronIterator(Electron const &block, bool begin);
  XDIAG_API ElectronIterator &operator++();
  XDIAG_API ProductState const &operator*() const;
  XDIAG_API bool operator!=(ElectronIterator const &rhs) const;

private:
  int64_t nsites_;
  mutable ProductState pstate_;
  basis::BasisElectronIterator it_;
};

} // namespace xdiag
