#pragma once

#include <variant>

#include <hydra/common.h>

#include <hydra/indexing/indexing_variants.h>

#include <hydra/blocks/spinhalf/terms/spinhalf_symmetric_exchange.h>
#include <hydra/blocks/spinhalf/terms/spinhalf_symmetric_scalar_chirality.h>
#include <hydra/blocks/spinhalf/terms/spinhalf_symmetric_ising.h>
#include <hydra/blocks/spinhalf/terms/spinhalf_unsymmetric_exchange.h>
#include <hydra/blocks/spinhalf/terms/spinhalf_unsymmetric_scalar_chirality.h>
#include <hydra/blocks/spinhalf/terms/spinhalf_unsymmetric_ising.h>
#include <hydra/blocks/spinhalf/terms/spinhalf_unsymmetric_spsm.h>
#include <hydra/blocks/spinhalf/terms/spinhalf_unsymmetric_sz.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Filler>
void spinhalf_ising(BondList const &bonds, Couplings const &couplings,
                    indexing::SpinhalfIndexing<bit_t> const &indexing,
                    Filler &&fill) {

  using namespace indexing;

  // Multiple dispatch using std::variant / std::visit
  std::visit(overloaded{[&](SpinhalfIndexingSz<bit_t> const &idxr) {
                          spinhalf_unsymmetric_ising<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](SpinhalfIndexingNoSz<bit_t> const &idxr) {
                          spinhalf_unsymmetric_ising<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](SpinhalfSymmetricIndexingSz<bit_t> const &idxr) {
                          spinhalf_symmetric_ising<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](SpinhalfSymmetricIndexingNoSz<bit_t> const &idxr) {
                          spinhalf_symmetric_ising<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](auto idxr) {
                          lila::Log.err("Invalid indexing encountered");
                        }},
             indexing);
}

template <typename bit_t, typename coeff_t, class Filler>
void spinhalf_exchange(BondList const &bonds, Couplings const &couplings,
                       indexing::SpinhalfIndexing<bit_t> const &indexing,
                       Filler &&fill) {

  using namespace indexing;

  // Multiple dispatch using std::variant / std::visit
  std::visit(overloaded{[&](SpinhalfIndexingSz<bit_t> const &idxr) {
                          spinhalf_unsymmetric_exchange<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](SpinhalfIndexingNoSz<bit_t> const &idxr) {
                          spinhalf_unsymmetric_exchange<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](SpinhalfSymmetricIndexingSz<bit_t> const &idxr) {
                          spinhalf_symmetric_exchange<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](SpinhalfSymmetricIndexingNoSz<bit_t> const &idxr) {
                          spinhalf_symmetric_exchange<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](auto idxr) {
                          lila::Log.err("Invalid indexing encountered");
                        }},
             indexing);
}

template <typename bit_t, typename coeff_t, class Filler>
void spinhalf_scalar_chirality(
    BondList const &bonds, Couplings const &couplings,
    indexing::SpinhalfIndexing<bit_t> const &indexing, Filler &&fill) {

  using namespace indexing;

  // Multiple dispatch using std::variant / std::visit
  std::visit(overloaded{[&](SpinhalfIndexingSz<bit_t> const &idxr) {
                          spinhalf_unsymmetric_scalar_chirality<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](SpinhalfIndexingNoSz<bit_t> const &idxr) {
                          spinhalf_unsymmetric_scalar_chirality<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](SpinhalfSymmetricIndexingSz<bit_t> const &idxr) {
                          spinhalf_symmetric_scalar_chirality<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](SpinhalfSymmetricIndexingNoSz<bit_t> const &idxr) {
                          spinhalf_symmetric_scalar_chirality<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](auto idxr) {
                          lila::Log.err("Invalid indexing encountered");
                        }},
             indexing);
}

template <typename bit_t, typename coeff_t, class Filler>
void spinhalf_sz(BondList const &bonds, Couplings const &couplings,
                 indexing::SpinhalfIndexing<bit_t> const &indexing,
                 Filler &&fill) {
  using namespace indexing;

  // Multiple dispatch using std::variant / std::visit
  std::visit(overloaded{[&](SpinhalfIndexingSz<bit_t> const &idxr) {
                          spinhalf_unsymmetric_sz<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](SpinhalfIndexingNoSz<bit_t> const &idxr) {
                          spinhalf_unsymmetric_sz<bit_t, coeff_t>(
                              bonds, couplings, idxr, fill);
                        },
                        [&](SpinhalfSymmetricIndexingSz<bit_t> const &idxr) {
                          // lila::Log.err("Sz (symmetric) not implemented");
                        },
                        [&](SpinhalfSymmetricIndexingNoSz<bit_t> const &idxr) {
                          // lila::Log.err("Sz (symmetric) not implemented");
                        },
                        [&](auto idxr) {
                          lila::Log.err("Invalid indexing encountered");
                        }},
             indexing);
}

template <typename bit_t, typename coeff_t, class Filler>
void spinhalf_spsm(BondList const &bonds, Couplings const &couplings,
                   indexing::SpinhalfIndexing<bit_t> const &indexing_in,
                   indexing::SpinhalfIndexing<bit_t> const &indexing_out,
                   Filler &&fill, std::string spsm) {

  using namespace indexing;

  // Multiple dispatch using std::variant / std::visit
  std::visit(
      overloaded{[&](SpinhalfIndexingSz<bit_t> const &idxr_in,
                     SpinhalfIndexingSz<bit_t> const &idxr_out) {
                   spinhalf_unsymmetric_spsm<bit_t, coeff_t>(
                       bonds, couplings, idxr_in, idxr_out, fill, spsm);
                 },
                 [&](SpinhalfIndexingNoSz<bit_t> const &idxr_in,
                     SpinhalfIndexingNoSz<bit_t> const &idxr_out) {
                   spinhalf_unsymmetric_spsm<bit_t, coeff_t>(
                       bonds, couplings, idxr_in, idxr_out, fill, spsm);
                 },
                 [&](SpinhalfSymmetricIndexingSz<bit_t> const &idxr_in,
                     SpinhalfSymmetricIndexingSz<bit_t> const &idxr_out) {
                   // lila::Log.err("S+/S- (symmetric) not implemented");
                 },
                 [&](SpinhalfSymmetricIndexingNoSz<bit_t> const &idxr_in,
                     SpinhalfSymmetricIndexingNoSz<bit_t> const &idxr_out) {
                   // lila::Log.err("S+/S- (symmetric) not implemented");
                 },
                 [&](auto idxr_in, auto idx_out) {
                   lila::Log.err("Invalid indexing encountered");
                 }},
      indexing_in, indexing_out);
}

} // namespace hydra::terms
