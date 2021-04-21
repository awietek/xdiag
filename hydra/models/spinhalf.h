#pragma once

#include <functional>

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>
#include <hydra/symmetries/representation.h>
#include <hydra/symmetries/spacegroup.h>

namespace hydra {

template <class bit_t> class SpinhalfIterator;

template <class bit_t, class SymmetryGroup = SpaceGroup<bit_t>>
class Spinhalf {
public:
  using iterator_t = SpinhalfIterator<bit_t>;

  Spinhalf() = default;
  explicit Spinhalf(int n_sites);
  Spinhalf(int n_sites, int sz);
  Spinhalf(int n_sites, SymmetryGroup symmetry_group, Representation irrep);
  Spinhalf(int n_sites, int sz, SymmetryGroup symmetry_group,
           Representation irrep);

  int n_sites() const { return n_sites_; }

  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }
  idx_t size() const { return size_; }

private:
  int n_sites_;
  bool sz_conserved_;
  int sz_;
  int n_upspins_;

  LinTable<bit_t> lintable_;

  bool symmetry_group_defined_;
  SymmetryGroup symmetry_group_;
  Representation irrep_;
  std::vector<bit_t> states_;
  std::vector<double> norms_;

  idx_t size_;
  iterator_t begin_;
  iterator_t end_;
};

// SpinhalfIterator
  template <class bit_t = std_bit_t> class SpinhalfIterator {
public:
  SpinhalfIterator() = default;

  template <class Advancer>
  SpinhalfIterator(bit_t state, idx_t index, Advancer &&advance)
      : state_(state), index_(index), advance_(advance) {}

  inline bool operator==(SpinhalfIterator<bit_t> const &rhs) const {
    return index_ == rhs.index_;
  }
  inline bool operator!=(SpinhalfIterator<bit_t> const &rhs) const {
    return !operator==(rhs);
  }
  inline SpinhalfIterator &operator++() {
    advance_(state_, index_);
    return *this;
  }
  inline bit_t operator*() const { return state_; }

private:
  bit_t state_;
  idx_t index_;
  std::function<void(bit_t &, idx_t &)> advance_;
};

} // namespace hydra
