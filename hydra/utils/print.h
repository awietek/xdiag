#pragma once

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <hydra/symmetries/permutation.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

#include <hydra/blocks/electron/electron.h>
#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/tj/tj.h>

#include <hydra/states/state.h>

#include <hydra/linalg/lanczos/tmatrix.h>
#include <lila/all.h>

namespace hydra::utils {

void PrintPretty(const char *identifier, double number);
void PrintPretty(const char *identifier, complex number);
  
void PrintPretty(const char *identifier, Bond const &bond);
void PrintPretty(const char *identifier, BondList const &bondlist);
void PrintPretty(const char *identifier, Couplings const &couplings);
void PrintPretty(const char *identifier, Permutation const &perm);
void PrintPretty(const char *identifier, PermutationGroup const &group);
void PrintPretty(const char *identifier, Representation const &irrep);

template <typename bit_t>
void PrintPretty(const char *identifier, Spinhalf<bit_t> const &block);
template <typename bit_t>
void PrintPretty(const char *identifier, tJ<bit_t> const &block);
template <typename bit_t>
void PrintPretty(const char *identifier, Electron<bit_t> const &block);

template <typename T>
void PrintPretty(const char *identifier, lila::Vector<T> const &vector) {
  lila::PrintPretty(identifier, vector);
}
template <typename T>
void PrintPretty(const char *identifier, lila::Matrix<T> const &matrix) {
  lila::PrintPretty(identifier, matrix);
}

void PrintPretty(const char *identifier, Tmatrix const &tmat);

template <typename coeff_t, class Block>
void PrintPretty(const char *identifier, State<coeff_t, Block> const &state) {
  printf("%s:\n", identifier);
  if constexpr (is_real<coeff_t>()) {
    printf("  real state\n");
  } else {
    printf("  cplx state\n");
  }
  printf("  size: %ld\n", state.size());
  printf("  norm: %.9e\n", lila::Norm(state.vector()));
  PrintPretty((std::string(identifier) + std::string(".block()")).c_str(),
              state.block());
}

} // namespace hydra::utils
