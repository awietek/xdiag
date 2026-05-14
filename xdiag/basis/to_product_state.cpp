#include "to_product_state.hpp"

#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>

namespace xdiag::basis {

template <typename enumeration_t>
ProductState to_product_state(enumeration_t const &enumeration, int64_t idx,
                              std::vector<std::string> const &dict) {
  int64_t nsites = enumeration.n();
  ProductState ps(nsites);
  auto const b = enumeration[idx];
  for (int64_t i = 0; i < nsites; ++i) {
    ps[i] = dict[bits::get(b, i)];
  }
  return ps;
}

#define INSTANTIATE_TO_PRODUCT_STATE(ENUMERATION_TYPE)                         \
  template ProductState to_product_state(ENUMERATION_TYPE const &, int64_t,    \
                                         std::vector<std::string> const &);

using namespace combinatorics;
using namespace bits;
INSTANTIATE_TO_PRODUCT_STATE(Subsets<uint32_t>);
INSTANTIATE_TO_PRODUCT_STATE(Subsets<uint64_t>);
INSTANTIATE_TO_PRODUCT_STATE(Combinations<uint32_t>);
INSTANTIATE_TO_PRODUCT_STATE(Combinations<uint64_t>);
INSTANTIATE_TO_PRODUCT_STATE(Combinations<BitsetStatic2>);
INSTANTIATE_TO_PRODUCT_STATE(Combinations<BitsetStatic4>);
INSTANTIATE_TO_PRODUCT_STATE(Combinations<BitsetStatic8>);
INSTANTIATE_TO_PRODUCT_STATE(Combinations<BitsetDynamic>);
INSTANTIATE_TO_PRODUCT_STATE(LinTable<uint32_t>);
INSTANTIATE_TO_PRODUCT_STATE(LinTable<uint64_t>);
} // namespace xdiag::basis
