#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {
template <class coeff_t, class Block>
coeff_t Inner(BondList const &bonds, Couplings const &couplings, Block &&block,
              lila::Vector<coeff_t> const &v) {
  auto Hv = lila::ZerosLike(v);
  Apply(bonds, couplings, block, v, block, Hv);
  return Dot(v, Hv);
}
} // namespace hydra
