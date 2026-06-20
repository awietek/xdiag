// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <utility>
#include <variant>

#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::utils {

// Helper type for visitor patterns
template <class... Ts> struct overload : Ts... {
  using Ts::operator()...;
};
template <class... Ts> overload(Ts...) -> overload<Ts...>;

template <class VariantA, class VariantB, class F>
decltype(auto) visit_same_type(VariantA &&a, VariantB &&b, F &&f,
                               std::string mismatch_message = "") try {
  if (a.index() != b.index()) {
    if (mismatch_message == "") {
      XDIAG_THROW("Type mismatch");
    } else {
      XDIAG_THROW(mismatch_message);
    }
  }
  return std::visit(
      [&](auto &&x) -> decltype(auto) {
        using T = std::decay_t<decltype(x)>;
        return std::forward<F>(f)(std::forward<decltype(x)>(x),
                                  std::get<T>(std::forward<VariantB>(b)));
      },
      std::forward<VariantA>(a));
}
XDIAG_CATCH

template <typename V1, typename V2> struct variant_union;

template <typename... Ts, typename... Us>
struct variant_union<std::variant<Ts...>, std::variant<Us...>> {
  using type = std::variant<Ts..., Us...>;
};

} // namespace xdiag::utils
