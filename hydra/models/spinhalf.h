#pragma once

#include <functional>

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>
#include <hydra/parameters/parameters.h>
#include <hydra/symmetries/charactertable.h>

namespace hydra {

  template <class bit_t> class SpinhalfIterator;

  
template <class bit_t = std_bit_t> class Spinhalf {
public:
  using iterator_t = SpinhalfIterator<bit_t>;

  Spinhalf() = default;
  Spinhalf(int n_sites);
  Spinhalf(int n_sites, Parameters p);
  Spinhalf(int n_sites, CharacterTable character_table, Parameters p);

  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }
  
  int n_sites() const { return n_sites_; }
  idx_t size() const { return size_; }

private:
  int n_sites_;
  Parameters parameters_;

  bool sz_conserved_;
  int sz_;
  LinTable<bit_t> lintable_;

  bool spinflip_symmetric_;
  double spinflip_sign_;

  bool space_group_defined_;
  std::string space_group_irrep_;
  CharacterTable character_table_;
  std::vector<complex> characters_;
  SpaceGroup space_group_;

  std::vector<bit_t> states_;
  std::vector<double> norms_;

  idx_t size_;
  iterator_t begin_;
  iterator_t end_;
};

// SpinhalfIterator
template <class bit_t> class SpinhalfIterator {
public:

  SpinhalfIterator() = default;

  template <class Advancer>
  SpinhalfIterator(bit_t state, idx_t index, Advancer && advance)
    : state_(state), index_(index), advance_(advance)
  {}

  inline bool operator==(SpinhalfIterator<bit_t> const &rhs) const {
    return index_ == rhs.index_;
  }
  inline bool operator!=(SpinhalfIterator<bit_t> const &rhs) const {
    return !operator==(rhs);
  }
  inline SpinhalfIterator &operator++() {
    state_ = advance_(state_);
    ++index_;
    return *this;
  }
  inline bit_t operator*() const { return state_; }

private:
  std::function<void(bit_t&, idx_t&)> advance_;
  bit_t state_;
  idx_t index_;
};

} // namespace hydra
