#pragma once

#include <memory>

#include <xdiag/basis/spinhalf/basis_spinhalf.hpp>
#include <xdiag/common.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

class SpinhalfIterator;

class Spinhalf {
public:
  using basis_t = basis::BasisSpinhalf;
  using iterator_t = SpinhalfIterator;

  Spinhalf() = default;
  Spinhalf(int64_t n_sites);
  Spinhalf(int64_t n_sites, int64_t n_up);
  Spinhalf(int64_t n_sites, PermutationGroup permutation_group,
           Representation irrep);
  Spinhalf(int64_t n_sites, int64_t n_up, PermutationGroup permutation_group,
           Representation irrep);

  // Constructors with sublattice coding
  Spinhalf(int64_t n_sites, PermutationGroup permutation_group,
           Representation irrep, int64_t n_sublat);
  Spinhalf(int64_t n_sites, int64_t n_up, PermutationGroup permutation_group,
           Representation irrep, int64_t n_sublat);

  int64_t n_sites() const;
  int64_t n_up() const;
  PermutationGroup permutation_group() const;
  Representation irrep() const;

  int64_t dim() const;
  int64_t size() const;
  iterator_t begin() const;
  iterator_t end() const;
  bool isreal(double precision = 1e-12) const;

  bool operator==(Spinhalf const &rhs) const;
  bool operator!=(Spinhalf const &rhs) const;

  basis_t const &basis() const;

private:
  int64_t n_sites_;
  int64_t n_up_;
  PermutationGroup permutation_group_;
  Representation irrep_;

  std::shared_ptr<basis_t> basis_;
  int64_t size_;
};

std::ostream &operator<<(std::ostream &out, Spinhalf const &block);
std::string to_string(Spinhalf const &block);

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
