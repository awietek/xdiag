// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <string_view>

#include <xdiag/bits/get_set.hpp>
#include <xdiag/states/product_state.hpp>

namespace xdiag::basis {

std::size_t create_basis_type_id();

class Basis {
public:
  virtual std::size_t type() const = 0;
  virtual std::string_view name() const = 0;
  virtual int64_t size() const = 0;
  virtual ProductState
  product_state(int64_t idx, std::vector<std::string> const &dict) const = 0;
  virtual ~Basis() = default;
};

template <typename Derived> class BasisType : public Basis {
public:
  static std::size_t static_type() {
    static const std::size_t id = create_basis_type_id();
    return id;
  }
  std::size_t type() const override { return static_type(); }
  std::string_view name() const override { return Derived::type_name; }

protected:
  ~BasisType() = default;
};

} // namespace xdiag::basis
