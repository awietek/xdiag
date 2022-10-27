#include "print.h"

#include <hydra/random/hashes.h>
#include <sstream>

namespace hydra::utils {

void PrintPretty(const char *identifier, std::string str) {
  printf("%s:\n", identifier);
  std::cout << str << "\n";
}

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

  if (bond.type_defined()) {
    printf("  type: %s\n", bond.type().c_str());
  } else {
    if (bond.is_real()) {
      PrintPretty("  matrix", bond.matrix_real());
    } else {
      PrintPretty("  matrix", bond.matrix());
    }
  }
  if (bond.coupling_defined()) {
    complex cpl = bond.coupling();

    if (std::abs(imag(cpl)) < 1e-12) {
      printf("  coupling: %.13e", real(cpl));
    } else {
      if (imag(cpl) > 0.) {
        printf("  coupling: %.17e + %.17eI\n", real(cpl), std::imag(cpl));
      } else {
        printf("  coupling: %.17e - %.17eI\n", real(cpl), -std::imag(cpl));
      }
    }
  } else {
    printf("  coupling_name: %s\n", bond.coupling_name().c_str());
  }

  printf("  sites: ");
  for (auto site : bond.sites())
    printf("%d ", site);
  printf("\n");
}

void PrintPretty(const char *identifier, BondList const &bondlist) {
  printf("%s:\n", identifier);
  for (auto bond : bondlist) {
    PrintPretty("bond", bond);
  }
}

void PrintPretty(const char *identifier, Permutation const &p) {
  printf("%s:\n  ", identifier);
  for (int i = 0; i < p.n_sites(); ++i) {
    printf("%d ", p[i]);
  }
  printf("\n");
  printf("  ID: 0x%lx\n", random::hash(p));
}

void PrintPretty(const char *identifier, PermutationGroup const &group) {
  printf("%s:\n", identifier);
  printf("  n_sites      : %d\n", group.n_sites());
  printf("  n_symmetries : %d\n", group.n_symmetries());
  printf("  ID           : 0x%lx\n", random::hash(group));
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
  printf("  ID        : 0x%lx\n", random::hash(irrep));
}

void PrintPretty(const char *identifier, Block const &block) {
  std::visit(
      overloaded{
          [identifier](Spinhalf const &block) {
            PrintPretty(identifier, block);
          },
          [identifier](tJ const &block) { PrintPretty(identifier, block); },
          [identifier](Electron const &block) {
            PrintPretty(identifier, block);
          },
      },
      block);
}

void PrintPretty(const char *identifier, Spinhalf const &block) {
  printf("%s:\n", identifier);

  printf("  n_sites  : %d\n", block.n_sites());
  if (block.sz_conserved()) {
    printf("  n_up     : %d\n", block.n_up());
  } else {
    printf("  n_up     : not conserved\n");
  }

  if (block.symmetric()) {
    printf("  group    : defined with ID 0x%lx\n",
           random::hash(block.permutation_group()));
    printf("  irrep    : defined with ID 0x%lx\n", random::hash(block.irrep()));
  }
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << block.size();
  printf("  dimension: %s\n", ss.str().c_str());
  printf("  ID       : 0x%lx\n", random::hash(block));
}

void PrintPretty(const char *identifier, tJ const &block) {
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
    printf("  group    : defined with ID 0x%lx\n",
           random::hash(block.permutation_group()));
    printf("  irrep    : defined with ID 0x%lx\n", random::hash(block.irrep()));
  }
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << block.size();
  printf("  dimension: %s\n", ss.str().c_str());
  printf("  ID       : 0x%lx\n", random::hash(block));
}

void PrintPretty(const char *identifier, Electron const &block) {
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
    printf("  group    : defined with ID 0x%lx\n",
           random::hash(block.permutation_group()));
    printf("  irrep    : defined with ID 0x%lx\n", random::hash(block.irrep()));
  }
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << block.size();
  printf("  dimension: %s\n", ss.str().c_str());
  printf("  ID       : 0x%lx\n", random::hash(block));
}

void PrintPretty(const char *identifier, RandomState const &rstate) {
  printf("%s:\n", identifier);
  printf("  RandomState, seed  : 0x%lx\n", rstate.seed());
}

void PrintPretty(const char *identifier, ProductState const &pstate) {
  printf("%s:\n", identifier);
  printf("  ProductState, n_sites: %d\n", pstate.n_sites());
  for (auto s : pstate) {
    printf("%s ", s.c_str());
  }
  printf("\n");
}

void PrintPretty(const char *identifier, ProductState const &pstate);

void PrintPretty(const char *identifier, Tmatrix const &tmat) {
  printf("%s:\n", identifier);
  PrintPretty("alphas", tmat.alphas());
  PrintPretty("betas", tmat.betas());
  PrintPretty("eigenvalues", tmat.eigenvalues());
}

void PrintPretty(const char *identifier, std::vector<float> const &vec) {
  PrintPretty(identifier, arma::Col<float>(vec));
}
void PrintPretty(const char *identifier, std::vector<double> const &vec) {
  PrintPretty(identifier, arma::Col<double>(vec));
}
void PrintPretty(const char *identifier, std::vector<scomplex> const &vec) {
  PrintPretty(identifier, arma::Col<scomplex>(vec));
}
void PrintPretty(const char *identifier, std::vector<complex> const &vec) {
  PrintPretty(identifier, arma::Col<complex>(vec));
}
void PrintPretty(const char *identifier, const arma::Mat<float> &mat) {
  printf("%s:\n", identifier);
  // for (uint32_t i = 0; i < mat.n_rows; ++i) {
  //   for (uint32_t j = 0; j < mat.n_cols; ++j)
  //     printf("%10.8g ", mat(i, j));
  //   printf("\n");
  // }
  // printf("\n");
  mat.brief_print();
}

void PrintPretty(const char *identifier, const arma::Mat<double> &mat) {
  printf("%s:\n", identifier);
  // for (uint32_t i = 0; i < mat.n_rows; ++i) {
  //   for (uint32_t j = 0; j < mat.n_cols; ++j)
  //     printf("%10.8g ", mat(i, j));
  //   printf("\n");
  // }
  // printf("\n");
  mat.brief_print();
}

void PrintPretty(const char *identifier,
                 const arma::Mat<std::complex<float>> &mat) {
  printf("%s:\n", identifier);
  //   for (uint32_t i = 0; i < mat.n_rows; ++i) {
  //     for (uint32_t j = 0; j < mat.n_cols; ++j)
  //       printf("%10.8g%-+8.8gj ", mat(i, j).real(), mat(i, j).imag());
  //     printf("\n");
  //   }
  //   printf("\n");
  mat.brief_print();
}

void PrintPretty(const char *identifier,
                 const arma::Mat<std::complex<double>> &mat) {
  printf("%s:\n", identifier);
  // for (uint32_t i = 0; i < mat.n_rows; ++i) {
  //   for (uint32_t j = 0; j < mat.n_cols; ++j)
  //     printf("%10.8g%-+8.8gj ", mat(i, j).real(), mat(i, j).imag());
  //   printf("\n");
  // }
  // printf("\n");
  mat.brief_print();
}

void PrintPretty(const char *identifier, const arma::Col<float> &vec) {
  printf("%s:\n", identifier);
  // for (uint32_t i = 0; i < vec.size(); ++i)
  //   printf("%10.8g ", vec(i));
  // printf("\n");
  vec.brief_print();
}

void PrintPretty(const char *identifier, const arma::Col<double> &vec) {
  printf("%s:\n", identifier);
  // for (uint32_t i = 0; i < vec.size(); ++i)
  //   printf("%10.8g ", vec(i));
  // printf("\n");
  vec.brief_print();
}

void PrintPretty(const char *identifier,
                 const arma::Col<std::complex<float>> &vec) {
  printf("%s:\n", identifier);
  // for (uint32_t i = 0; i < vec.size(); ++i)
  //   printf("%10.8g%-+8.8gj ", vec(i).real(), vec(i).imag());
  // printf("\n");
  vec.brief_print();
}

void PrintPretty(const char *identifier,
                 const arma::Col<std::complex<double>> &vec) {
  printf("%s:\n", identifier);
  // for (uint32_t i = 0; i < vec.size(); ++i)
  //   printf("%10.8g%-+8.8gj ", vec(i).real(), vec(i).imag());
  // printf("\n");
  vec.brief_print();
}

} // namespace hydra::utils
