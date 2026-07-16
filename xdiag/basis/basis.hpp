// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <memory>
#include <string_view>
#include <utility>

#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/states/product_state.hpp>

namespace xdiag::basis {

// Sentinel returned by index lookups when a configuration is not present in the
// (possibly rank-local) basis.
inline constexpr int64_t invalid_index = -1;

std::size_t create_basis_type_id();

// Type-erased forward iterator over a basis, yielding ProductStates. This is
// the only abstraction needed to iterate a type-erased Basis lazily: one step
// (advance) and one dereference (product_state). The block iterators drive it.
class BasisIterator {
public:
  virtual ~BasisIterator() = default;
  virtual void advance() = 0;
  virtual ProductState product_state() const = 0;
};

class Basis {
public:
  virtual std::size_t type() const = 0;
  virtual std::string_view name() const = 0;
  virtual int64_t size() const = 0;
  virtual int64_t nsites() const = 0;

  // next three functions are for distributed blocks which override this default
  // implementation
  virtual int64_t dim() const { return size(); }
  virtual int64_t size_max() const { return size(); }
  virtual int64_t size_min() const { return size(); }
  virtual std::unique_ptr<BasisIterator> product_state_iterator() const = 0;
  virtual int64_t index(ProductState const &pstate) const = 0;
  virtual ~Basis() = default;
};

// Per-site local-state index of a basis configuration, used when converting a
// raw configuration to a ProductState. A single bit_t (spin-1/2 / boson /
// single fermion species backends) gives bits::get(config, i) directly. A
// spinful electron basis dereferences to a (ups, dns) pair, whose local state
// is ups_i + 2 * dns_i (0 empty, 1 up, 2 dn, 3 up&dn).
template <typename bit_t>
inline int64_t local_state(bit_t const &config, int64_t i) {
  return bits::get(config, i);
}
template <typename bit_t>
inline int64_t local_state(std::pair<bit_t, bit_t> const &config, int64_t i) {
  return bits::get(config.first, i) + 2 * bits::get(config.second, i);
}

// Inverse of local_state: reconstruct a raw configuration from a ProductState,
// used by the Basis::index overrides. bits::zero sizes the storage per backend.
template <typename bit_t>
inline bit_t config_from_pstate(ProductState const &pstate, int64_t nsites) {
  bit_t config = bits::zero<bit_t>(nsites);
  for (int64_t i = 0; i < nsites; ++i) {
    bits::set(config, i, pstate[i]);
  }
  return config;
}

// (ups, dns) product backends: local state 0 empty, 1 up, 2 dn, 3 up&dn.
template <typename bit_t>
inline std::pair<bit_t, bit_t>
pair_from_pstate(ProductState const &pstate, int64_t nsites) {
  bit_t ups = bits::zero<bit_t>(nsites);
  bit_t dns = bits::zero<bit_t>(nsites);
  for (int64_t i = 0; i < nsites; ++i) {
    int64_t l = pstate[i];
    bits::set(ups, i, l & 1);
    bits::set(dns, i, (l >> 1) & 1);
  }
  return {ups, dns};
}

// Concrete BasisIterator wrapping a basis' own linear iterator. The underlying
// iterator dereferences to the raw configuration; the conversion to a
// ProductState happens here at the leaf, where the configuration type is known.
// local_state extracts the per-site local-state index, so this works for any
// local dimension and for the (ups, dns) electron product configuration.
template <typename Derived> class BasisIteratorImpl : public BasisIterator {
public:
  BasisIteratorImpl(typename Derived::iterator_t it, int64_t nsites)
      : it_(it), nsites_(nsites) {}

  void advance() override { ++it_; }

  ProductState product_state() const override {
    auto config = *it_;
    ProductState ps(nsites_);
    for (int64_t i = 0; i < nsites_; ++i) {
      ps[i] = local_state(config, i);
    }
    return ps;
  }

private:
  typename Derived::iterator_t it_;
  int64_t nsites_;
};

template <typename Derived> class BasisType : public Basis {
public:
  static std::size_t static_type() {
    static const std::size_t id = create_basis_type_id();
    return id;
  }
  std::size_t type() const override { return static_type(); }
  std::string_view name() const override { return Derived::type_name; }

  std::unique_ptr<BasisIterator> product_state_iterator() const override {
    Derived const *derived = static_cast<Derived const *>(this);
    return std::make_unique<BasisIteratorImpl<Derived>>(derived->begin(),
                                                        derived->nsites());
  }

protected:
  ~BasisType() = default;
};

} // namespace xdiag::basis
