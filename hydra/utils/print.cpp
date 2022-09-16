#include "print.h"

#include <hydra/random/hashes.h>
#include <sstream>

namespace hydra::utils {

void PrintPretty(const char *identifier, int number) {
  printf("%s:\n", identifier);
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << number;
  printf("%s\n", ss.str().c_str());
}

void PrintPretty(const char *identifier, uint32_t number) {
  printf("%s:\n", identifier);
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << number;
  printf("%s\n", ss.str().c_str());
}

void PrintPretty(const char *identifier, uint64_t number) {
  printf("%s:\n", identifier);
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << number;
  printf("%s\n", ss.str().c_str());
}

void PrintPretty(const char *identifier, int64_t number) {
  printf("%s:\n", identifier);
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << number;
  printf("%s\n", ss.str().c_str());
}

void PrintPretty(const char *identifier, double number) {
  printf("%s:\n", identifier);
  printf("%.17e\n", number);
}

void PrintPretty(const char *identifier, complex number) {
  printf("%s:\n", identifier);
  if (std::imag(number) > 0.) {
    printf("%.17e + %.17eI\n", std::real(number), std::imag(number));
  } else {
    printf("%.17e - %.17eI\n", std::real(number), -std::imag(number));
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
  printf("  size      : %ld\n", (long)irrep.size());
  printf("  characters:");
  for (auto c : irrep.characters()) {
    if (std::imag(c) > 0.) {
      printf("(%.9f + %.9fI) ", std::real(c), std::imag(c));
    } else if (std::imag(c) == 0.) {
      printf("(%.9f + %.9fI) ", std::real(c), 0.);
    } else {
      printf("(%.9f - %.9fI) ", std::real(c), -std::imag(c));
    }
  }
  printf("\n");
  printf("  ID        : 0x%08x\n", random::hash(irrep));
}

template <typename bit_t>
void PrintPretty(const char *identifier, Spinhalf<bit_t> const &block) {
  printf("%s:\n", identifier);

  printf("  n_sites  : %d\n", block.n_sites());
  if (block.sz_conserved()) {
    printf("  n_up     : %d\n", block.n_up());
  } else {
    printf("  n_up     : not conserved\n");
  }

  if (block.symmetric()) {
    printf("  group    : defined with ID 0x%08x\n",
           random::hash(block.permutation_group()));
    printf("  irrep    : defined with ID 0x%08x\n",
           random::hash(block.irrep()));
  }
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << block.size();
  printf("  dimension: %s\n", ss.str().c_str());
  printf("  ID       : 0x%08x\n", random::hash(block));
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

  printf("  n_sites  : %d\n", block.n_sites());
  if (block.sz_conserved() && block.charge_conserved()) {
    printf("  n_up     : %d\n", block.n_up());
    printf("  n_dn     : %d\n", block.n_dn());

  } else {
    printf("  n_up     : not conserved\n");
    printf("  n_dn     : not conserved\n");
  }

  if (block.symmetric()) {
    printf("  group    : defined with ID 0x%08x\n",
           random::hash(block.permutation_group()));
    printf("  irrep    : defined with ID 0x%08x\n",
           random::hash(block.irrep()));
  }
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << block.size();
  printf("  dimension: %s\n", ss.str().c_str());
  printf("  ID       : 0x%08x\n", random::hash(block));
}
template void PrintPretty(const char *identifier, tJ<uint16_t> const &block);
template void PrintPretty(const char *identifier, tJ<uint32_t> const &block);
template void PrintPretty(const char *identifier, tJ<uint64_t> const &block);

template <typename bit_t>
void PrintPretty(const char *identifier, Electron<bit_t> const &block) {
  printf("%s:\n", identifier);

  printf("  n_sites  : %d\n", block.n_sites());
  if (block.sz_conserved() && block.charge_conserved()) {
    printf("  n_up     : %d\n", block.n_up());
    printf("  n_dn     : %d\n", block.n_dn());

  } else {
    printf("  n_up     : not conserved\n");
    printf("  n_dn     : not conserved\n");
  }

  if (block.symmetric()) {
    printf("  group    : defined with ID 0x%08x\n",
           random::hash(block.permutation_group()));
    printf("  irrep    : defined with ID 0x%08x\n",
           random::hash(block.irrep()));
  }
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << block.size();
  printf("  dimension: %s\n", ss.str().c_str());
  printf("  ID       : 0x%08x\n", random::hash(block));
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

void PrintPretty(const char *identifier, const arma::Mat<float> &mat) {
  printf("%s:\n", identifier);
  for (uint32_t i = 0; i < mat.n_rows; ++i) {
    for (uint32_t j = 0; j < mat.n_cols; ++j)
      printf("%10.8g ", mat(i, j));
    printf("\n");
  }
  printf("\n");
}

void PrintPretty(const char *identifier, const arma::Mat<double> &mat) {
  printf("%s:\n", identifier);
  for (uint32_t i = 0; i < mat.n_rows; ++i) {
    for (uint32_t j = 0; j < mat.n_cols; ++j)
      printf("%10.8g ", mat(i, j));
    printf("\n");
  }
  printf("\n");
}

void PrintPretty(const char *identifier,
                 const arma::Mat<std::complex<float>> &mat) {
  printf("%s:\n", identifier);
  for (uint32_t i = 0; i < mat.n_rows; ++i) {
    for (uint32_t j = 0; j < mat.n_cols; ++j)
      printf("%10.8g%-+8.8gj ", mat(i, j).real(), mat(i, j).imag());
    printf("\n");
  }
  printf("\n");
}

void PrintPretty(const char *identifier,
                 const arma::Mat<std::complex<double>> &mat) {
  printf("%s:\n", identifier);
  for (uint32_t i = 0; i < mat.n_rows; ++i) {
    for (uint32_t j = 0; j < mat.n_cols; ++j)
      printf("%10.8g%-+8.8gj ", mat(i, j).real(), mat(i, j).imag());
    printf("\n");
  }
  printf("\n");
}

void PrintPretty(const char *identifier, const arma::Col<float> &vec) {
  printf("%s:\n", identifier);
  for (uint32_t i = 0; i < vec.size(); ++i)
    printf("%10.8g ", vec(i));
  printf("\n");
}

void PrintPretty(const char *identifier, const arma::Col<double> &vec) {
  printf("%s:\n", identifier);
  for (uint32_t i = 0; i < vec.size(); ++i)
    printf("%10.8g ", vec(i));
  printf("\n");
}

void PrintPretty(const char *identifier,
                 const arma::Col<std::complex<float>> &vec) {
  printf("%s:\n", identifier);
  for (uint32_t i = 0; i < vec.size(); ++i)
    printf("%10.8g%-+8.8gj ", vec(i).real(), vec(i).imag());
  printf("\n");
}

void PrintPretty(const char *identifier,
                 const arma::Col<std::complex<double>> &vec) {
  printf("%s:\n", identifier);
  for (uint32_t i = 0; i < vec.size(); ++i)
    printf("%10.8g%-+8.8gj ", vec(i).real(), vec(i).imag());
  printf("\n");
}

} // namespace hydra::utils
