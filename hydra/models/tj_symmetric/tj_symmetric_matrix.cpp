#include "tj_symmetric_matrix.h"

#include <hydra/utils/bitops.h>

#include <hydra/models/tj_symmetric/terms/tj_symmetric_hopping.h>
#include <hydra/models/tj_symmetric/terms/tj_symmetric_exchange.h>
#include <hydra/models/tj_symmetric/terms/tj_symmetric_ising.h>

namespace hydra {

template <class bit_t, class SymmetryGroup>
lila::Matrix<complex>
MatrixCplx(BondList const &bonds, Couplings const &couplings,
            tJSymmetric<bit_t, SymmetryGroup> const &block_in,
            tJSymmetric<bit_t, SymmetryGroup> const &block_out) {
  assert(block_in == block_out); // only temporary

  idx_t dim = block_in.size();
  auto mat = lila::Zeros<complex>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, complex val) {
    mat(idx_out, idx_in) += val;
  };

  tj::do_hopping_symmetric<bit_t, complex>(bonds, couplings, block_in, fill);
  tj::do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  tj::do_exchange_symmetric<bit_t, complex>(bonds, couplings, block_in, fill);

  return mat;
}

template lila::Matrix<complex> MatrixCplx<uint16, SpaceGroup<uint16>>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetric<uint16, SpaceGroup<uint16>> const &block_in,
    tJSymmetric<uint16, SpaceGroup<uint16>> const &block_out);
template lila::Matrix<complex> MatrixCplx<uint32, SpaceGroup<uint32>>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetric<uint32, SpaceGroup<uint32>> const &block_in,
    tJSymmetric<uint32, SpaceGroup<uint32>> const &block_out);
template lila::Matrix<complex> MatrixCplx<uint64, SpaceGroup<uint64>>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetric<uint64, SpaceGroup<uint64>> const &block_in,
    tJSymmetric<uint64, SpaceGroup<uint64>> const &block_out);

} // namespace hydra
