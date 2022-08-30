#include "print.h"

#include <hydra/random/hashes.h>

namespace hydra::utils {

void PrintPretty(const char *identifier, double number) {
  printf("%s:\n", identifier);
  printf("%.17e\n", number);
}

void PrintPretty(const char *identifier, complex number) {
  printf("%s:\n", identifier);
  if (std::imag(number) > 0.) {
    printf("%.17e + %.17e I\n", std::real(number), std::imag(number));
  } else {
    printf("%.17e - %.17e I\n", std::real(number), -std::imag(number));
  }
}

void PrintPretty(const char *identifier, Bond const &bond) {
  printf("%s:\n", identifier);
  printf("%s %s ", bond.type().c_str(), bond.coupling().c_str());
  for (auto site : bond.sites())
    printf("%d ", site);
  printf("\n");
}

void PrintPretty(const char *identifier, BondList const &bondlist) {
  printf("%s:\n", identifier);
  for (auto bond : bondlist) {
    printf("%s %s ", bond.type().c_str(), bond.coupling().c_str());
    for (auto site : bond.sites())
      printf("%d ", site);
    printf("\n");
  }
}

void PrintPretty(const char *identifier, Couplings const &couplings) {
  printf("%s:\n", identifier);
  for (auto coupling : couplings) {
    printf("%s %f %f\n", coupling.first.c_str(), std::real(coupling.second),
           std::imag(coupling.second));
  }
}

void PrintPretty(const char *identifier, Permutation const &p) {
  printf("%s:\n  ", identifier);
  for (int i = 0; i < p.n_sites(); ++i) {
    printf("%d ", p[i]);
  }
  printf("\n");
  printf("  ID: 0x%08x\n", random::hash(p));
}

void PrintPretty(const char *identifier, PermutationGroup const &group) {
  printf("%s:\n", identifier);
  printf("  n_sites      : %d\n", group.n_sites());
  printf("  n_symmetries : %d\n", group.n_symmetries());
  printf("  ID           : 0x%08x\n", random::hash(group));
}

void PrintPretty(const char *identifier, Representation const &irrep) {
  printf("%s:\n", identifier);
  printf("  size : %ld\n", irrep.size());
  printf("  ID   : 0x%08x\n", random::hash(irrep));
}

template <typename bit_t>
void PrintPretty(const char *identifier, Spinhalf<bit_t> const &block) {
  printf("%s:\n", identifier);

  printf("  n_sites: %d\n", block.n_sites());
  if (block.sz_conserved()) {
    printf("  n_up   : %d\n", block.n_up());
  } else {
    printf("  n_up   : not conserved\n");
  }

  if (block.symmetric()) {
    printf("  group  : defined with ID 0x%08x\n",
           random::hash(block.permutation_group()));
    printf("  irrep  : defined with ID 0x%08x\n", random::hash(block.irrep()));
  }

  printf("  ID     : 0x%08x\n", random::hash(block));
}
template void PrintPretty(const char *identifier,
                          Spinhalf<uint16_t> const &block);
template void PrintPretty(const char *identifier,
                          Spinhalf<uint32_t> const &block);
template void PrintPretty(const char *identifier,
                          Spinhalf<uint64_t> const &block);

template <typename bit_t>
void PrintPretty(const char *identifier, tJ<bit_t> const &block) {
  printf("%s:\n", identifier);

  printf("  n_sites: %d\n", block.n_sites());
  if (block.sz_conserved() && block.charge_conserved()) {
    printf("  n_up   : %d\n", block.n_up());
    printf("  n_dn   : %d\n", block.n_dn());

  } else {
    printf("  n_up   : not conserved\n");
    printf("  n_dn   : not conserved\n");
  }

  if (block.symmetric()) {
    printf("  group  : defined with ID 0x%08x\n",
           random::hash(block.permutation_group()));
    printf("  irrep  : defined with ID 0x%08x\n", random::hash(block.irrep()));
  }

  printf("  ID     : 0x%08x\n", random::hash(block));
}
template void PrintPretty(const char *identifier, tJ<uint16_t> const &block);
template void PrintPretty(const char *identifier, tJ<uint32_t> const &block);
template void PrintPretty(const char *identifier, tJ<uint64_t> const &block);

template <typename bit_t>
void PrintPretty(const char *identifier, Electron<bit_t> const &block) {
  printf("%s:\n", identifier);

  printf("  n_sites: %d\n", block.n_sites());
  if (block.sz_conserved() && block.charge_conserved()) {
    printf("  n_up   : %d\n", block.n_up());
    printf("  n_dn   : %d\n", block.n_dn());

  } else {
    printf("  n_up   : not conserved\n");
    printf("  n_dn   : not conserved\n");
  }

  if (block.symmetric()) {
    printf("  group  : defined with ID 0x%08x\n",
           random::hash(block.permutation_group()));
    printf("  irrep  : defined with ID 0x%08x\n", random::hash(block.irrep()));
  }

  printf("  ID     : 0x%08x\n", random::hash(block));
}
template void PrintPretty(const char *identifier,
                          Electron<uint16_t> const &block);
template void PrintPretty(const char *identifier,
                          Electron<uint32_t> const &block);
template void PrintPretty(const char *identifier,
                          Electron<uint64_t> const &block);

template <typename bit_t>
void PrintPretty(const char *identifier, tJ<bit_t> const &block);
template <typename bit_t>
void PrintPretty(const char *identifier, Electron<bit_t> const &block);

void PrintPretty(const char *identifier, Tmatrix const &tmat) {
  printf("%s:\n", identifier);
  PrintPretty("alphas", tmat.alphas());
  PrintPretty("betas", tmat.betas());
  PrintPretty("eigenvalues", tmat.eigenvalues());
}

} // namespace hydra::utils
