#include "electron_symmetric_apply.h"

// #include <hydra/models/electron_symmetric/terms/electron_symmetric_exchange.h>
#include <hydra/models/electron_symmetric/terms/electron_symmetric_hopping.h>
// #include <hydra/models/electron_symmetric/terms/electron_symmetric_ising.h>
#include <hydra/models/electron_symmetric/terms/electron_symmetric_u.h>

#include <hydra/models/utils/model_utils.h>

namespace hydra {

template <class bit_t, class GroupAction>
void Apply(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetric<bit_t, GroupAction> const &block_in,
           lila::Vector<double> const &vec_in,
           ElectronSymmetric<bit_t, GroupAction> const &block_out,
           lila::Vector<double> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  utils::check_symmetric_operator_real(bonds, couplings, block_in.irrep(),
                                       block_out.irrep(),
                                       "apply real ElectronSymmetric operator");

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, double val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  electronterms::do_U_symmetric(couplings, block_in, fill);
  electronterms::do_hopping_symmetric<bit_t, double>(bonds, couplings, block_in,
                                                     fill);
  // electronterms::do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  // electronterms::do_exchange_symmetric<bit_t, double>(bonds, couplings,
  //                                                     block_in, fill);
}

template void Apply<uint16, PermutationGroupLookup<uint16>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_in,
    lila::Vector<double> const &vec_in,
    ElectronSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_out,
    lila::Vector<double> &vec_out);
template void Apply<uint32, PermutationGroupLookup<uint32>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_in,
    lila::Vector<double> const &vec_in,
    ElectronSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_out,
    lila::Vector<double> &vec_out);
template void Apply<uint64, PermutationGroupLookup<uint64>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_in,
    lila::Vector<double> const &vec_in,
    ElectronSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_out,
    lila::Vector<double> &vec_out);

template <class bit_t, class GroupAction>
void Apply(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetric<bit_t, GroupAction> const &block_in,
           lila::Vector<complex> const &vec_in,
           ElectronSymmetric<bit_t, GroupAction> const &block_out,
           lila::Vector<complex> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, complex val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  electronterms::do_U_symmetric(couplings, block_in, fill);
  electronterms::do_hopping_symmetric<bit_t, complex>(bonds, couplings,
                                                      block_in, fill);
  // electronterms::do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  // electronterms::do_exchange_symmetric<bit_t, complex>(bonds, couplings,
  //                                                      block_in, fill);
}

template void Apply<uint16, PermutationGroupLookup<uint16>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_in,
    lila::Vector<complex> const &vec_in,
    ElectronSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint32, PermutationGroupLookup<uint32>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_in,
    lila::Vector<complex> const &vec_in,
    ElectronSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint64, PermutationGroupLookup<uint64>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_in,
    lila::Vector<complex> const &vec_in,
    ElectronSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_out,
    lila::Vector<complex> &vec_out);

} // namespace hydra
