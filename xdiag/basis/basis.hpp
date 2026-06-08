// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <memory>
#include <string_view>

#include <xdiag/bits/get_set.hpp>
#include <xdiag/states/product_state.hpp>

namespace xdiag::basis {

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
  // Iterator positioned at the first basis element. Advancing it walks the
  // basis linearly; nothing else is needed since end is detected by index.
  virtual std::unique_ptr<BasisIterator> product_state_iterator() const = 0;
  virtual ~Basis() = default;
};

// Concrete BasisIterator wrapping a basis' own linear iterator. The underlying
// iterator dereferences to the raw bit configuration; the conversion to a
// ProductState happens here at the leaf, where the bit type is known. bits::get
// extracts the per-site local-state index (a single bit for spin-1/2 backends,
// an nbits field for BitArray backends), so this works for any local dimension.
template <typename Derived> class BasisIteratorImpl : public BasisIterator {
public:
  BasisIteratorImpl(typename Derived::iterator_t it, int64_t nsites)
      : it_(it), nsites_(nsites) {}

  void advance() override { ++it_; }

  ProductState product_state() const override {
    typename Derived::bit_t config = *it_;
    ProductState ps(nsites_);
    for (int64_t i = 0; i < nsites_; ++i) {
      ps[i] = bits::get(config, i);
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
