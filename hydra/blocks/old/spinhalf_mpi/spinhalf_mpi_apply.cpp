#include "spinhalf_mpi_apply.h"

#include <hydra/bitops/bitops.h>
#include <hydra/utils/logger.h>
#include <hydra/combinatorics/combinations.h>

#include <hydra/operators/operator_qns.h>

#include <hydra/blocks/spinhalf_mpi/terms/spinhalf_mpi_terms.h>
#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <class bit_t, class coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           SpinhalfMPI<bit_t> const &block_in,
           lila::Vector<coeff_t> const &vec_in,
           SpinhalfMPI<bit_t> const &block_out,
           lila::Vector<coeff_t> &vec_out) {

  int n_up_out = utils::spinhalf_nup(bonds, couplings, block_in);
  if (n_up_out != block_out.n_up())
    Log.err("Incompatible n_up in Apply: {} != {}", n_up_out,
	    block_out.n_up());

  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  utils::check_operator_works_with<coeff_t>(bonds, couplings, "spinhalf_apply");

  lila::Zeros(vec_out);

  auto const &indexing_in = block_in.indexing();
  auto const &indexing_out = block_out.indexing();

  if (block_in == block_out) {

    auto tis = rightnow_mpi();
    terms::spinhalf_mpi_ising(bonds, couplings, indexing_in, vec_in, vec_out);
    timing_mpi(tis, rightnow_mpi(), " (ising)", 2);

    auto tex = rightnow_mpi();
    terms::spinhalf_mpi_exchange(bonds, couplings, indexing_in, vec_in,
                                 vec_out);
    timing_mpi(tex, rightnow_mpi(), " (exchange)", 2);

    auto tsz = rightnow_mpi();
    terms::spinhalf_mpi_sz(bonds, couplings, indexing_in, vec_in, vec_out);
    timing_mpi(tsz, rightnow_mpi(), " (sz)", 2);
  }
  auto tspsm = rightnow_mpi();
  terms::spinhalf_mpi_spsm(bonds, couplings, indexing_in, vec_in, indexing_out,
                           vec_out, "S+");
  terms::spinhalf_mpi_spsm(bonds, couplings, indexing_in, vec_in, indexing_out,
                           vec_out, "S-");
  timing_mpi(tspsm, rightnow_mpi(), " (spsm)", 2);
}

template void Apply<uint16_t>(BondList const &bonds, Couplings const &couplings,
                              SpinhalfMPI<uint16_t> const &block_in,
                              lila::Vector<double> const &vec_in,
                              SpinhalfMPI<uint16_t> const &block_out,
                              lila::Vector<double> &vec_out);
template void Apply<uint32_t>(BondList const &bonds, Couplings const &couplings,
                              SpinhalfMPI<uint32_t> const &block_in,
                              lila::Vector<double> const &vec_in,
                              SpinhalfMPI<uint32_t> const &block_out,
                              lila::Vector<double> &vec_out);
template void Apply<uint64_t>(BondList const &bonds, Couplings const &couplings,
                              SpinhalfMPI<uint64_t> const &block_in,
                              lila::Vector<double> const &vec_in,
                              SpinhalfMPI<uint64_t> const &block_out,
                              lila::Vector<double> &vec_out);

template void Apply<uint16_t, complex>(BondList const &bonds,
                                       Couplings const &couplings,
                                       SpinhalfMPI<uint16_t> const &block_in,
                                       lila::Vector<complex> const &vec_in,
                                       SpinhalfMPI<uint16_t> const &block_out,
                                       lila::Vector<complex> &vec_out);
template void Apply<uint32_t, complex>(BondList const &bonds,
                                       Couplings const &couplings,
                                       SpinhalfMPI<uint32_t> const &block_in,
                                       lila::Vector<complex> const &vec_in,
                                       SpinhalfMPI<uint32_t> const &block_out,
                                       lila::Vector<complex> &vec_out);
template void Apply<uint64_t, complex>(BondList const &bonds,
                                       Couplings const &couplings,
                                       SpinhalfMPI<uint64_t> const &block_in,
                                       lila::Vector<complex> const &vec_in,
                                       SpinhalfMPI<uint64_t> const &block_out,
                                       lila::Vector<complex> &vec_out);

} // namespace hydra
