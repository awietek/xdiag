#include "electron_symmetric_simple_apply.h"

#include <hydra/models/electron_symmetric_simple/terms/electron_symmetric_simple_exchange.h>
#include <hydra/models/electron_symmetric_simple/terms/electron_symmetric_simple_hopping.h>
#include <hydra/models/electron_symmetric_simple/terms/electron_symmetric_simple_ising.h>
#include <hydra/models/electron_symmetric_simple/terms/electron_symmetric_simple_u.h>

#include <hydra/models/utils/model_utils.h>

namespace hydra {

template <class bit_t, class GroupAction>
void Apply(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetricSimple<bit_t, GroupAction> const &block_in,
           lila::Vector<double> const &vec_in,
           ElectronSymmetricSimple<bit_t, GroupAction> const &block_out,
           lila::Vector<double> &vec_out) {
  using namespace terms::electron_symmetric_simple;

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  utils::check_symmetric_operator_real(
      bonds, couplings, block_in.irrep(), block_out.irrep(),
      "apply real ElectronSymmetricSimple operator");

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, double val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  do_U_symmetric(couplings, block_in, fill);
  do_hopping_symmetric<bit_t, double>(bonds, couplings, block_in, fill);
  do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  do_exchange_symmetric<bit_t, double>(bonds, couplings, block_in, fill);
}

template void Apply<uint16, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetricSimple<uint16, PermutationGroupAction> const &block_in,
    lila::Vector<double> const &vec_in,
    ElectronSymmetricSimple<uint16, PermutationGroupAction> const &block_out,
    lila::Vector<double> &vec_out);
template void Apply<uint32, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetricSimple<uint32, PermutationGroupAction> const &block_in,
    lila::Vector<double> const &vec_in,
    ElectronSymmetricSimple<uint32, PermutationGroupAction> const &block_out,
    lila::Vector<double> &vec_out);
template void Apply<uint64, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetricSimple<uint64, PermutationGroupAction> const &block_in,
    lila::Vector<double> const &vec_in,
    ElectronSymmetricSimple<uint64, PermutationGroupAction> const &block_out,
    lila::Vector<double> &vec_out);

template <class bit_t, class GroupAction>
void Apply(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetricSimple<bit_t, GroupAction> const &block_in,
           lila::Vector<complex> const &vec_in,
           ElectronSymmetricSimple<bit_t, GroupAction> const &block_out,
           lila::Vector<complex> &vec_out) {
  using namespace terms::electron_symmetric_simple;

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, complex val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  do_U_symmetric(couplings, block_in, fill);
  do_hopping_symmetric<bit_t, complex>(bonds, couplings, block_in, fill);
  do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  do_exchange_symmetric<bit_t, complex>(bonds, couplings, block_in, fill);
}

template void Apply<uint16, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetricSimple<uint16, PermutationGroupAction> const &block_in,
    lila::Vector<complex> const &vec_in,
    ElectronSymmetricSimple<uint16, PermutationGroupAction> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint32, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetricSimple<uint32, PermutationGroupAction> const &block_in,
    lila::Vector<complex> const &vec_in,
    ElectronSymmetricSimple<uint32, PermutationGroupAction> const &block_out,
    lila::Vector<complex> &vec_out);
template void Apply<uint64, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetricSimple<uint64, PermutationGroupAction> const &block_in,
    lila::Vector<complex> const &vec_in,
    ElectronSymmetricSimple<uint64, PermutationGroupAction> const &block_out,
    lila::Vector<complex> &vec_out);

} // namespace hydra
