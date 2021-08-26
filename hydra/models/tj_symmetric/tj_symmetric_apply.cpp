#include "tj_symmetric_apply.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/tj_symmetric/terms/tj_symmetric_hopping.h>
#include <hydra/models/tj_symmetric/terms/tj_symmetric_exchange.h>
#include <hydra/models/tj_symmetric/terms/tj_symmetric_ising.h>

namespace hydra {

template <class bit_t, class GroupAction>
void Apply(BondList const &bonds, Couplings const &couplings,
           tJSymmetric<bit_t, GroupAction> const &block_in,
           lila::Vector<double> const &vec_in,
           tJSymmetric<bit_t, GroupAction> const &block_out,
           lila::Vector<double> &vec_out) {
  if (is_complex(bonds))
    lila::Log.err("Cannot apply to real vector with complex bonds!");
  if (is_complex(couplings))
    lila::Log.err("Cannot apply to real vector with complex couplings!");
  if (is_complex(block_in.irrep()) || is_complex(block_out.irrep()))
    lila::Log.err("Cannot apply to real vector with complex representation!");

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, double val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  tj::do_hopping_symmetric<bit_t, double>(bonds, couplings, block_in, fill);
  tj::do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  tj::do_exchange_symmetric<bit_t, double>(bonds, couplings, block_in, fill);
}

template void Apply<uint16, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetric<uint16, PermutationGroupAction> const &block_in,
    lila::Vector<double> const &vec_in,
    tJSymmetric<uint16, PermutationGroupAction> const &block_out,
    lila::Vector<double> &vec_out);
template void Apply<uint32, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetric<uint32, PermutationGroupAction> const &block_in,
    lila::Vector<double> const &vec_in,
    tJSymmetric<uint32, PermutationGroupAction> const &block_out,
    lila::Vector<double> &vec_out);
template void Apply<uint64, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetric<uint64, PermutationGroupAction> const &block_in,
    lila::Vector<double> const &vec_in,
    tJSymmetric<uint64, PermutationGroupAction> const &block_out,
    lila::Vector<double> &vec_out);


template <class bit_t, class GroupAction>
void Apply(BondList const &bonds, Couplings const &couplings,
           tJSymmetric<bit_t, GroupAction> const &block_in,
           lila::Vector<complex> const &vec_in,
           tJSymmetric<bit_t, GroupAction> const &block_out,
           lila::Vector<complex> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, complex val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  tj::do_hopping_symmetric<bit_t, complex>(bonds, couplings, block_in, fill);
  tj::do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  tj::do_exchange_symmetric<bit_t, complex>(bonds, couplings, block_in, fill);
}

template void Apply<uint16, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetric<uint16, PermutationGroupAction> const &block_in,
    lila::Vector<complex> const &vec_in,
    tJSymmetric<uint16, PermutationGroupAction> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint32, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetric<uint32, PermutationGroupAction> const &block_in,
    lila::Vector<complex> const &vec_in,
    tJSymmetric<uint32, PermutationGroupAction> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint64, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    tJSymmetric<uint64, PermutationGroupAction> const &block_in,
    lila::Vector<complex> const &vec_in,
    tJSymmetric<uint64, PermutationGroupAction> const &block_out,
    lila::Vector<complex> &vec_out);

} // namespace hydra
