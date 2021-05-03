#include "electron_symmetric_matrix.h"

#include <hydra/utils/bitops.h>

#include <hydra/models/electron/electron_utils.h>
#include <hydra/models/electron/terms/electron_symmetric_u.h>
#include <hydra/models/electron/terms/electron_symmetric_hopping.h>

namespace hydra {

template <class bit_t, class SymmetryGroup>
lila::Matrix<complex>
matrix_cplx(BondList const &bonds, Couplings const &couplings,
            ElectronSymmetric<bit_t, SymmetryGroup> const &block_in,
            ElectronSymmetric<bit_t, SymmetryGroup> const &block_out) {
  assert(block_in == block_out); // only temporary

  idx_t dim = block_in.size();
  auto mat = lila::Zeros<complex>(dim, dim);

  electron::do_U_symmetric(couplings, block_in, [&mat](idx_t idx, double val) {
    mat(idx, idx) += val;
  });
  
  electron::do_hopping_symmetric<bit_t, complex>(
      bonds, couplings, block_in,
      [&mat](idx_t idx_out, idx_t idx_in, complex val) {
        mat(idx_out, idx_in) += val;
      });

  return mat;
}

template lila::Matrix<complex> matrix_cplx<uint16, SpaceGroup<uint16>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint16, SpaceGroup<uint16>> const &block_in,
    ElectronSymmetric<uint16, SpaceGroup<uint16>> const &block_out);
template lila::Matrix<complex> matrix_cplx<uint32, SpaceGroup<uint32>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint32, SpaceGroup<uint32>> const &block_in,
    ElectronSymmetric<uint32, SpaceGroup<uint32>> const &block_out);
template lila::Matrix<complex> matrix_cplx<uint64, SpaceGroup<uint64>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint64, SpaceGroup<uint64>> const &block_in,
    ElectronSymmetric<uint64, SpaceGroup<uint64>> const &block_out);

} // namespace hydra
