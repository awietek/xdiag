#include "tj_symmetric_simple_apply.h"

#include <hydra/blocks/tj_symmetric_simple/terms/tj_symmetric_simple_exchange.h>
#include <hydra/blocks/tj_symmetric_simple/terms/tj_symmetric_simple_hopping.h>
#include <hydra/blocks/tj_symmetric_simple/terms/tj_symmetric_simple_ising.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <class bit_t, class GroupAction>
void Apply(BondList const &bonds, Couplings const &couplings,
           tJSymmetricSimple<bit_t, GroupAction> const &block_in,
           lila::Vector<double> const &vec_in,
           tJSymmetricSimple<bit_t, GroupAction> const &block_out,
           lila::Vector<double> &vec_out) {
  using namespace terms::tj_symmetric_simple;

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  utils::check_symmetric_operator_real(bonds, couplings, block_in.irrep(),
                                       block_out.irrep(),
                                       "apply real tJSymmetricSimple operator");

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, double val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  do_hopping_symmetric<bit_t, double>(bonds, couplings, block_in, fill);
  do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  do_exchange_symmetric<bit_t, double>(bonds, couplings, block_in, fill);
}

template void Apply<uint16_t, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetricSimple<uint16_t, PermutationGroupAction> const &block_in,
    lila::Vector<double> const &vec_in,
    tJSymmetricSimple<uint16_t, PermutationGroupAction> const &block_out,
    lila::Vector<double> &vec_out);
template void Apply<uint32_t, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetricSimple<uint32_t, PermutationGroupAction> const &block_in,
    lila::Vector<double> const &vec_in,
    tJSymmetricSimple<uint32_t, PermutationGroupAction> const &block_out,
    lila::Vector<double> &vec_out);
template void Apply<uint64_t, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetricSimple<uint64_t, PermutationGroupAction> const &block_in,
    lila::Vector<double> const &vec_in,
    tJSymmetricSimple<uint64_t, PermutationGroupAction> const &block_out,
    lila::Vector<double> &vec_out);

template <class bit_t, class GroupAction>
void Apply(BondList const &bonds, Couplings const &couplings,
           tJSymmetricSimple<bit_t, GroupAction> const &block_in,
           lila::Vector<complex> const &vec_in,
           tJSymmetricSimple<bit_t, GroupAction> const &block_out,
           lila::Vector<complex> &vec_out) {
  using namespace terms::tj_symmetric_simple;

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, complex val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  do_hopping_symmetric<bit_t, complex>(bonds, couplings, block_in, fill);
  do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  do_exchange_symmetric<bit_t, complex>(bonds, couplings, block_in, fill);
}

template void Apply<uint16_t, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetricSimple<uint16_t, PermutationGroupAction> const &block_in,
    lila::Vector<complex> const &vec_in,
    tJSymmetricSimple<uint16_t, PermutationGroupAction> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint32_t, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetricSimple<uint32_t, PermutationGroupAction> const &block_in,
    lila::Vector<complex> const &vec_in,
    tJSymmetricSimple<uint32_t, PermutationGroupAction> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint64_t, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetricSimple<uint64_t, PermutationGroupAction> const &block_in,
    lila::Vector<complex> const &vec_in,
    tJSymmetricSimple<uint64_t, PermutationGroupAction> const &block_out,
    lila::Vector<complex> &vec_out);

} // namespace hydra
