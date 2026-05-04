#include "basis_symmetric.hpp"

#include <xdiag/basis/to_product_state.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/symmetries/action/site_permutation.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::basis {

template <typename enumeration_t>
BasisSymmetric<enumeration_t>::BasisSymmetric(enumeration_t const &enumeration,
                                              Representation const &irrep) try
    : enumeration_(enumeration), irrep_(irrep), table_(enumeration, irrep) {}
XDIAG_CATCH

template <typename enumeration_t>
int64_t BasisSymmetric<enumeration_t>::size() const {
  return table_.size();
}

template <typename enumeration_t>
int64_t BasisSymmetric<enumeration_t>::nsites() const {
  return enumeration_.n();
}

template <typename enumeration_t>
Representation const &BasisSymmetric<enumeration_t>::irrep() const {
  return irrep_;
}

template <typename enumeration_t>
typename BasisSymmetric<enumeration_t>::iterator_t
BasisSymmetric<enumeration_t>::begin() const {
  return table_.begin();
}

template <typename enumeration_t>
typename BasisSymmetric<enumeration_t>::iterator_t
BasisSymmetric<enumeration_t>::end() const {
  return table_.end();
}

template <typename enumeration_t>
ProductState BasisSymmetric<enumeration_t>::product_state(
    int64_t idx, std::vector<std::string> const &dict) const {
  return to_product_state(enumeration_, idx, dict);
}

template <typename enumeration_t>
bool BasisSymmetric<enumeration_t>::operator==(
    BasisSymmetric<enumeration_t> const &rhs) const {
  return (irrep_ == rhs.irrep_) && (enumeration_ == rhs.enumeration_);
}

template <typename enumeration_t>
bool BasisSymmetric<enumeration_t>::operator!=(
    BasisSymmetric<enumeration_t> const &rhs) const {
  return !operator==(rhs);
}

using namespace combinatorics;
using namespace bits;
template class BasisSymmetric<Subsets<uint32_t>>;
template class BasisSymmetric<Subsets<uint64_t>>;
template class BasisSymmetric<Combinations<uint32_t>>;
template class BasisSymmetric<Combinations<uint64_t>>;
template class BasisSymmetric<Combinations<BitsetStatic2>>;
template class BasisSymmetric<Combinations<BitsetStatic4>>;
template class BasisSymmetric<Combinations<BitsetStatic8>>;
template class BasisSymmetric<Combinations<BitsetDynamic>>;
template class BasisSymmetric<LinTable<uint32_t>>;
template class BasisSymmetric<LinTable<uint64_t>>;

} // namespace xdiag::basis
