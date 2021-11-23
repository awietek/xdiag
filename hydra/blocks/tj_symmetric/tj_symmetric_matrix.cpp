#include "tj_symmetric_matrix.h"

#include <hydra/blocks/tj_symmetric/terms/tj_symmetric_exchange.h>
#include <hydra/blocks/tj_symmetric/terms/tj_symmetric_hopping.h>
#include <hydra/blocks/tj_symmetric/terms/tj_symmetric_ising.h>
#include <hydra/blocks/tj_symmetric/tj_symmetric.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <class bit_t, typename coeff_t>
lila::Matrix<coeff_t> MatrixGen(BondList const &bonds,
                                Couplings const &couplings,
                                tJSymmetric<bit_t> const &block_in,
                                tJSymmetric<bit_t> const &block_out) {
  using namespace terms::tj_symmetric;

  assert(block_in == block_out); // only temporary

  utils::check_operator_works_with<coeff_t>(bonds, couplings, block_in.irrep(),
                                            block_out.irrep(),
                                            "tj_symmetric_apply");
  idx_t dim = block_in.size();
  auto mat = lila::Zeros<coeff_t>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, coeff_t val) {
    mat(idx_out, idx_in) += val;
  };

  auto const &indexing_in = block_in.indexing();
  // auto const &indexing_out = block_out.indexing();

  do_hopping_symmetric<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  do_ising_symmetric<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  do_exchange_symmetric<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);

  return mat;
}

template lila::Matrix<double>
MatrixGen<uint16_t, double>(BondList const &, Couplings const &,
                            tJSymmetric<uint16_t> const &,
                            tJSymmetric<uint16_t> const &);
template lila::Matrix<double>
MatrixGen<uint32_t, double>(BondList const &, Couplings const &,
                            tJSymmetric<uint32_t> const &,
                            tJSymmetric<uint32_t> const &);
template lila::Matrix<double>
MatrixGen<uint64_t, double>(BondList const &, Couplings const &,
                            tJSymmetric<uint64_t> const &,
                            tJSymmetric<uint64_t> const &);

template lila::Matrix<complex>
MatrixGen<uint16_t, complex>(BondList const &, Couplings const &,
                             tJSymmetric<uint16_t> const &,
                             tJSymmetric<uint16_t> const &);
template lila::Matrix<complex>
MatrixGen<uint32_t, complex>(BondList const &, Couplings const &,
                             tJSymmetric<uint32_t> const &,
                             tJSymmetric<uint32_t> const &);
template lila::Matrix<complex>
MatrixGen<uint64_t, complex>(BondList const &, Couplings const &,
                             tJSymmetric<uint64_t> const &,
                             tJSymmetric<uint64_t> const &);

} // namespace hydra
