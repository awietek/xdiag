#include "spinhalf_symmetric_apply.h"

#include <hydra/models/spinhalf_symmetric/terms/spinhalf_symmetric_exchange.h>
#include <hydra/models/spinhalf_symmetric/terms/spinhalf_symmetric_ising.h>

#include <hydra/models/utils/model_utils.h>

namespace hydra {

template <class bit_t, class GroupAction>
void Apply(BondList const &bonds, Couplings const &couplings,
           SpinhalfSymmetric<bit_t, GroupAction> const &block_in,
           lila::Vector<double> const &vec_in,
           SpinhalfSymmetric<bit_t, GroupAction> const &block_out,
           lila::Vector<double> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  utils::check_symmetric_operator_real(bonds, couplings, block_in.irrep(),
                                       block_out.irrep(),
                                       "apply real SpinhalfSymmetric operator");

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, double val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  spinhalfterms::do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  spinhalfterms::do_exchange_symmetric<bit_t, double>(bonds, couplings,
                                                      block_in, fill);
}

template void Apply<uint16, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint16, PermutationGroupAction> const &block_in,
    lila::Vector<double> const &vec_in,
    SpinhalfSymmetric<uint16, PermutationGroupAction> const &block_out,
    lila::Vector<double> &vec_out);
template void Apply<uint32, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint32, PermutationGroupAction> const &block_in,
    lila::Vector<double> const &vec_in,
    SpinhalfSymmetric<uint32, PermutationGroupAction> const &block_out,
    lila::Vector<double> &vec_out);
template void Apply<uint64, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint64, PermutationGroupAction> const &block_in,
    lila::Vector<double> const &vec_in,
    SpinhalfSymmetric<uint64, PermutationGroupAction> const &block_out,
    lila::Vector<double> &vec_out);

template void Apply<uint16, PermutationGroupLookup<uint16>>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_in,
    lila::Vector<double> const &vec_in,
    SpinhalfSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_out,
    lila::Vector<double> &vec_out);
template void Apply<uint32, PermutationGroupLookup<uint32>>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_in,
    lila::Vector<double> const &vec_in,
    SpinhalfSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_out,
    lila::Vector<double> &vec_out);
template void Apply<uint64, PermutationGroupLookup<uint64>>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_in,
    lila::Vector<double> const &vec_in,
    SpinhalfSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_out,
    lila::Vector<double> &vec_out);

template <class bit_t, class GroupAction>
void Apply(BondList const &bonds, Couplings const &couplings,
           SpinhalfSymmetric<bit_t, GroupAction> const &block_in,
           lila::Vector<complex> const &vec_in,
           SpinhalfSymmetric<bit_t, GroupAction> const &block_out,
           lila::Vector<complex> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, complex val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  spinhalfterms::do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  spinhalfterms::do_exchange_symmetric<bit_t, complex>(bonds, couplings,
                                                       block_in, fill);
}

template void Apply<uint16, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint16, PermutationGroupAction> const &block_in,
    lila::Vector<complex> const &vec_in,
    SpinhalfSymmetric<uint16, PermutationGroupAction> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint32, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint32, PermutationGroupAction> const &block_in,
    lila::Vector<complex> const &vec_in,
    SpinhalfSymmetric<uint32, PermutationGroupAction> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint64, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint64, PermutationGroupAction> const &block_in,
    lila::Vector<complex> const &vec_in,
    SpinhalfSymmetric<uint64, PermutationGroupAction> const &block_out,
    lila::Vector<complex> &vec_out);

template void Apply<uint16, PermutationGroupLookup<uint16>>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_in,
    lila::Vector<complex> const &vec_in,
    SpinhalfSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint32, PermutationGroupLookup<uint32>>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_in,
    lila::Vector<complex> const &vec_in,
    SpinhalfSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint64, PermutationGroupLookup<uint64>>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_in,
    lila::Vector<complex> const &vec_in,
    SpinhalfSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_out,
    lila::Vector<complex> &vec_out);

} // namespace hydra
