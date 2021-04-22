#pragma once

#include <functional>

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>
#include <hydra/symmetries/fermi_sign.h>
#include <hydra/symmetries/representation.h>
#include <hydra/symmetries/spacegroup.h>

namespace hydra {

template <class bit_t> class ElectronIterator;

template <class bit_t, class SymmetryGroup = SpaceGroup<bit_t>> class Electron {
public:
  using iterator_t = ElectronIterator<bit_t>;

  Electron() = default;
  Electron(int n_sites, int charge, int sz);
  Electron(int n_sites, int charge, int sz, SymmetryGroup symmetry_group,
           Representation irrep);

  int n_sites() const { return n_sites_; }

  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }
  idx_t size() const { return size_; }

private:
  int n_sites_;

  bool charge_conserved_;
  int charge_;
  bool sz_conserved_;
  int sz_;
  int nup_;
  int ndn_;

  LinTable<bit_t> lintable_ups_;
  LinTable<bit_t> lintable_dns_;

  bool symmetry_group_defined_;
  SymmetryGroup symmetry_group_;
  Representation irrep_;
  std::vector<bit_t> upspins_;
  std::vector<bit_t> dnspins_;
  std::vector<double> norms_;

  idx_t size_;
  iterator_t begin_;
  iterator_t end_;
};

// ElectronIterator
template <class bit_t = std_bit_t> class ElectronIterator {
public:
  ElectronIterator() = default;

  template <class Advancer>
  ElectronIterator(bit_t ups, bit_t dns, idx_t index, Advancer &&advance)
      : ups_(ups), dns_(dns), index_(index), advance_(advance) {}

  inline bool operator==(ElectronIterator<bit_t> const &rhs) const {
    return index_ == rhs.index_;
  }
  inline bool operator!=(ElectronIterator<bit_t> const &rhs) const {
    return !operator==(rhs);
  }
  inline ElectronIterator &operator++() {
    advance_(ups_, dns_, index_);
    return *this;
  }
  inline std::tuple<bit_t, bit_t> operator*() const { return {ups_, dns_}; }
  inline bit_t ups() const { return ups_; }
  inline bit_t dns() const { return dns_; }

private:
  bit_t ups_;
  bit_t dns_;
  idx_t index_;
  std::function<void(bit_t &, bit_t &, idx_t &)> advance_;
};

} // namespace hydra
