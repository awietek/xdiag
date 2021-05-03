#include "electron_symmetric_apply.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/electron/terms/electron_symmetric_hopping.h>
#include <hydra/models/electron/terms/electron_symmetric_u.h>

namespace hydra {

template <class bit_t, class SymmetryGroup>
void apply(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetric<bit_t, SymmetryGroup> const &block_in,
           lila::Vector<complex> const &vec_in,
           ElectronSymmetric<bit_t, SymmetryGroup> const &block_out,
           lila::Vector<complex> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);
  electron::do_U_symmetric(couplings, block_in,
                           [&vec_out, &vec_in](idx_t idx, double val) {
                             vec_out(idx) += val * vec_in(idx);
                           });
  
  electron::do_hopping_symmetric<bit_t, complex>(
      bonds, couplings, block_in,
      [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, complex val) {
        vec_out(idx_out) += val * vec_in(idx_in);
      });
}

template void apply<uint16, SpaceGroup<uint16>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint16, SpaceGroup<uint16>> const &block_in,
    lila::Vector<complex> const &vec_in,
    ElectronSymmetric<uint16, SpaceGroup<uint16>> const &block_out,
    lila::Vector<complex> &vec_out);
template void apply<uint32, SpaceGroup<uint32>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint32, SpaceGroup<uint32>> const &block_in,
    lila::Vector<complex> const &vec_in,
    ElectronSymmetric<uint32, SpaceGroup<uint32>> const &block_out,
    lila::Vector<complex> &vec_out);
template void apply<uint64, SpaceGroup<uint64>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint64, SpaceGroup<uint64>> const &block_in,
    lila::Vector<complex> const &vec_in,
    ElectronSymmetric<uint64, SpaceGroup<uint64>> const &block_out,
    lila::Vector<complex> &vec_out);

} // namespace hydra
