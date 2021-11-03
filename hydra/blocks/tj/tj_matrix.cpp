#include "tj_matrix.h"

#include <hydra/blocks/tj/terms/tj_exchange.h>
#include <hydra/blocks/tj/terms/tj_hopping.h>
#include <hydra/blocks/tj/terms/tj_ising.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <class bit_t>
lila::Matrix<double>
MatrixReal(BondList const &bonds, Couplings const &couplings,
           tJ<bit_t> const &block_in, tJ<bit_t> const &block_out) {
  using namespace terms::tj;

  assert(block_in == block_out); // only temporary
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  utils::check_operator_real(bonds, couplings, "construct real tJ matrix");

  auto mat = lila::Zeros<double>(dim_out, dim_in);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, double val) {
    mat(idx_out, idx_in) += val;
  };

  do_hopping<bit_t, double>(bonds, couplings, block_in, fill);
  do_ising<bit_t>(bonds, couplings, block_in, fill);
  do_exchange<bit_t>(bonds, couplings, block_in, fill);
  return mat;
}

template <class bit_t>
lila::Matrix<complex>
MatrixCplx(BondList const &bonds, Couplings const &couplings,
           tJ<bit_t> const &block_in, tJ<bit_t> const &block_out) {
  using namespace terms::tj;

  assert(block_in == block_out);
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  auto mat = lila::Zeros<complex>(dim_out, dim_in);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, complex val) {
    mat(idx_out, idx_in) += val;
  };

  do_hopping<bit_t, complex>(bonds, couplings, block_in, fill);
  do_ising<bit_t>(bonds, couplings, block_in, fill);
  do_exchange<bit_t>(bonds, couplings, block_in, fill);
  return mat;
}

template lila::Matrix<double>
MatrixReal<uint16_t>(BondList const &bonds, Couplings const &couplings,
                     tJ<uint16_t> const &block_in,
                     tJ<uint16_t> const &block_out);
template lila::Matrix<double>
MatrixReal<uint32_t>(BondList const &bonds, Couplings const &couplings,
                     tJ<uint32_t> const &block_in,
                     tJ<uint32_t> const &block_out);
template lila::Matrix<double>
MatrixReal<uint64_t>(BondList const &bonds, Couplings const &couplings,
                     tJ<uint64_t> const &block_in,
                     tJ<uint64_t> const &block_out);

template lila::Matrix<complex>
MatrixCplx<uint16_t>(BondList const &bonds, Couplings const &couplings,
                     tJ<uint16_t> const &block_in,
                     tJ<uint16_t> const &block_out);
template lila::Matrix<complex>
MatrixCplx<uint32_t>(BondList const &bonds, Couplings const &couplings,
                     tJ<uint32_t> const &block_in,
                     tJ<uint32_t> const &block_out);
template lila::Matrix<complex>
MatrixCplx<uint64_t>(BondList const &bonds, Couplings const &couplings,
                     tJ<uint64_t> const &block_in,
                     tJ<uint64_t> const &block_out);

} // namespace hydra
