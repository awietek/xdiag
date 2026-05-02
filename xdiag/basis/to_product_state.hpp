#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <xdiag/states/product_state.hpp>

namespace xdiag::basis {

template <typename enumeration_t>
ProductState to_product_state(enumeration_t const &enumeration, int64_t idx,
                              std::vector<std::string> const &dict);

} // namespace xdiag::basis
