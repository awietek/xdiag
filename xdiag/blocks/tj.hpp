#pragma once
#include <xdiag/basis/tj/basis_tj.hpp>
#include <xdiag/common.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

class tJIterator;

class tJ {
public:
  using basis_t = basis::BasistJ;
  using iterator_t = tJIterator;

  tJ() = default;
  tJ(int64_t n_sites, int64_t nup, int64_t ndn);
  tJ(int64_t n_sites, int64_t nup, int64_t ndn,
     PermutationGroup permutation_group, Representation irrep);

  int64_t n_sites() const;
  int64_t n_up() const;
  int64_t n_dn() const;
  PermutationGroup const &permutation_group() const;
  Representation const &irrep() const;

  int64_t dim() const;
  int64_t size() const;
  iterator_t begin() const;
  iterator_t end() const;
  int64_t index(ProductState const& pstate) const;
  bool isreal(double precision = 1e-12) const;

  bool operator==(tJ const &rhs) const;
  bool operator!=(tJ const &rhs) const;

  basis_t const &basis() const;

private:
  int64_t n_sites_;
  int64_t n_up_;
  int64_t n_dn_;
  PermutationGroup permutation_group_;
  Representation irrep_;
  std::shared_ptr<basis_t> basis_;
  int64_t size_;
};

std::ostream &operator<<(std::ostream &out, tJ const &block);
std::string to_string(tJ const &block);

class tJIterator {
public:
  tJIterator(tJ const &block, bool begin);
  tJIterator &operator++();
  ProductState const &operator*() const;
  bool operator!=(tJIterator const &rhs) const;

private:
  int64_t n_sites_;
  mutable ProductState pstate_;
  basis::BasistJIterator it_;
};

} // namespace xdiag
