#include "electron_symmetric_apply.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/electron_symmetric/terms/electron_symmetric_hopping.h>
#include <hydra/models/electron_symmetric/terms/electron_symmetric_u.h>

namespace hydra {

template <class bit_t, class SymmetryGroup>
void Apply(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetric<bit_t, SymmetryGroup> const &block_in,
           lila::Vector<complex> const &vec_in,
           ElectronSymmetric<bit_t, SymmetryGroup> const &block_out,
           lila::Vector<complex> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, complex val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  electron::do_U_symmetric(couplings, block_in, fill);
  electron::do_hopping_symmetric<bit_t, complex>(bonds, couplings, block_in,
                                                 fill);
}

template void Apply<uint16, SpaceGroup<uint16>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint16, SpaceGroup<uint16>> const &block_in,
    lila::Vector<complex> const &vec_in,
    ElectronSymmetric<uint16, SpaceGroup<uint16>> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint32, SpaceGroup<uint32>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint32, SpaceGroup<uint32>> const &block_in,
    lila::Vector<complex> const &vec_in,
    ElectronSymmetric<uint32, SpaceGroup<uint32>> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint64, SpaceGroup<uint64>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint64, SpaceGroup<uint64>> const &block_in,
    lila::Vector<complex> const &vec_in,
    ElectronSymmetric<uint64, SpaceGroup<uint64>> const &block_out,
    lila::Vector<complex> &vec_out);

} // namespace hydra
