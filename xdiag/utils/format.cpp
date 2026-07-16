// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

// Compiles fmt's out-of-line core exactly once (compiled mode), instead of
// re-instantiating it header-only in every translation unit. All other sources
// include <extern/fmt/format.hpp> without FMT_HEADER_ONLY, getting declarations.
//
// The explicit instantiations below mirror fmt's own src/format.cc: they provide
// the char / std::locale specializations that fmt otherwise instantiates
// per-TU in header-only mode.
#include <extern/fmt/format-inl.hpp>

FMT_BEGIN_NAMESPACE
namespace detail {

template FMT_API auto dragonbox::to_decimal(float x) noexcept
    -> dragonbox::decimal_fp<float>;
template FMT_API auto dragonbox::to_decimal(double x) noexcept
    -> dragonbox::decimal_fp<double>;

#ifndef FMT_STATIC_THOUSANDS_SEPARATOR
template FMT_API locale_ref::locale_ref(const std::locale &loc);
template FMT_API auto locale_ref::get<std::locale>() const -> std::locale;
#endif

template FMT_API auto thousands_sep_impl(locale_ref)
    -> thousands_sep_result<char>;
template FMT_API auto decimal_point_impl(locale_ref) -> char;

template FMT_API void buffer<char>::append(const char *, const char *);

template FMT_API void vformat_to(buffer<char> &buf, string_view fmt,
                                 typename vformat_args<>::type args,
                                 locale_ref loc);

} // namespace detail
FMT_END_NAMESPACE
