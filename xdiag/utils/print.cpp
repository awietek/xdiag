#include "print.hpp"

#include <cinttypes>
#include <sstream>

#include <xdiag/random/hash.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::utils {

void print_pretty(const char *identifier, std::string str) {
  printf("%s:\n", identifier);
  std::cout << str << "\n";
}

void print_pretty(const char *identifier, int number) {
  printf("%s:\n", identifier);
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << number;
  printf("%s\n", ss.str().c_str());
}

void print_pretty(const char *identifier, uint32_t number) {
  printf("%s:\n", identifier);
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << number;
  printf("%s\n", ss.str().c_str());
}

void print_pretty(const char *identifier, uint64_t number) {
  printf("%s:\n", identifier);
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << number;
  printf("%s\n", ss.str().c_str());
}

void print_pretty(const char *identifier, int64_t number) {
  printf("%s:\n", identifier);
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << number;
  printf("%s\n", ss.str().c_str());
}

void print_pretty(const char *identifier, double number) {
  printf("%s:\n", identifier);
  printf("%.17e\n", number);
}

void print_pretty(const char *identifier, complex number) {
  printf("%s:\n", identifier);
  if (std::imag(number) > 0.) {
    printf("%.17e + %.17eI\n", std::real(number), std::imag(number));
  } else {
    printf("%.17e - %.17eI\n", std::real(number), -std::imag(number));
  }
}

void print_pretty(const char *identifier, Bond const &bond) {
  printf("%s:\n", identifier);
  std::cout << bond;
}

void print_pretty(const char *identifier, BondList const &bondlist) {
  printf("%s:\n", identifier);
  printf("Interactions:\n");
  for (auto bond : bondlist) {
    print_pretty(" bond", bond);
  }
  if (bondlist.couplings().size() > 0) {
    printf("Couplings:\n");
    for (auto name : bondlist.couplings()) {
      auto cpl = bondlist[name];
      std::visit([](auto &&arg) { std::cout << arg << "\n"; }, cpl);
    }
  }
}

void print_pretty(const char *identifier, Permutation const &p) {
  printf("%s:\n  ", identifier);
  for (int i = 0; i < p.size(); ++i) {
    printf("%" PRId64 " ", p[i]);
  }
  printf("\n");
  printf("  ID: 0x%lx\n", (unsigned long)random::hash(p));
}

void print_pretty(const char *identifier, PermutationGroup const &group) {
  printf("%s:\n", identifier);
  printf("  n_sites      : %" PRId64 "\n", group.n_sites());
  printf("  n_symmetries : %" PRId64 "\n", group.n_symmetries());
  printf("  ID           : 0x%lx\n", (unsigned long)random::hash(group));
}

void print_pretty(const char *identifier, Representation const &irrep) {
  printf("%s:\n", identifier);
  printf("  size      : %" PRId64 "\n", irrep.size());
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
  printf("  ID        : 0x%lx\n", (unsigned long)random::hash(irrep));
}

void print_pretty(const char *identifier, U1 const &g) {
  (void)g;
  printf("%s:\n", identifier);
  printf("  U(1) group\n");
}

void print_pretty(const char *identifier, QNum const &qn) {
  printf("%s:\n", identifier);
  std::visit(overload{
                 [](U1 g, int i) {
                   (void)g;
                   printf("  QNum, U(1), n=%d\n", i);
                 },
                 [](PermutationGroup const &g, Representation const &i) {
                   printf("  QNum\n");
                   printf("  PermutationGroup\n");
                   printf("    n_sites      : %" PRId64 "\n", g.n_sites());
                   printf("    n_symmetries : %" PRId64 "\n", g.n_symmetries());
                   printf("    ID           : 0x%lx\n\n",
                          (unsigned long)random::hash(g));
                   printf("  Representation\n");
                   printf("    ID        : 0x%lx\n",
                          (unsigned long)random::hash(i));
                 },
                 [](auto const &g, auto const &i) {
                   (void)g;
                   (void)i;
                   throw std::runtime_error(
                       "Error printing QNum: incompatible group or irrep");
                 },
             },
             qn.group(), qn.irrep());
}
void print_pretty(const char *identifier, QN const &qn) {
  (void)identifier;
  (void)qn;
}

void print_pretty(const char *identifier, block_variant_t const &block) {
  printf("llllll\n");
  std::visit(
      overload{
          [identifier](Spinhalf const &block) {
            print_pretty(identifier, block);
          },
          [identifier](tJ const &block) { print_pretty(identifier, block); },
          [identifier](Electron const &block) {
            print_pretty(identifier, block);
          },
#ifdef XDIAG_USE_MPI
          [identifier](tJDistributed const &block) {
            print_pretty(identifier, block);
          },
#endif
      },
      block);
}

void print_pretty(const char *identifier, Spinhalf const &block) {
  printf("%s:\n", identifier);

  printf("  n_sites  : %" PRId64 "\n", block.n_sites());
  if (block.sz_conserved()) {
    printf("  n_up     : %" PRId64 "\n", block.n_up());
  } else {
    printf("  n_up     : not conserved\n");
  }

  if (block.symmetric()) {
    printf("  group    : defined with ID 0x%lx\n",
           (unsigned long)random::hash(block.permutation_group()));
    printf("  irrep    : defined with ID 0x%lx\n",
           (unsigned long)random::hash(block.irrep()));
  }
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << block.size();
  printf("  dimension: %s\n", ss.str().c_str());
  printf("  ID       : 0x%lx\n", (unsigned long)random::hash(block));
}

void print_pretty(const char *identifier, tJ const &block) {
  printf("%s:\n", identifier);

  printf("  n_sites  : %" PRId64 "\n", block.n_sites());
  if (block.sz_conserved() && block.charge_conserved()) {
    printf("  n_up     : %" PRId64 "\n", block.n_up());
    printf("  n_dn     : %" PRId64 "\n", block.n_dn());

  } else {
    printf("  n_up     : not conserved\n");
    printf("  n_dn     : not conserved\n");
  }

  if (block.symmetric()) {
    printf("  group    : defined with ID 0x%lx\n",
           (unsigned long)random::hash(block.permutation_group()));
    printf("  irrep    : defined with ID 0x%lx\n",
           (unsigned long)random::hash(block.irrep()));
  }
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << block.size();
  printf("  dimension: %s\n", ss.str().c_str());
  printf("  ID       : 0x%lx\n", (unsigned long)random::hash(block));
}

#ifdef XDIAG_USE_MPI
void print_pretty(const char *identifier, tJDistributed const &block) {

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (rank == 0) {
    printf("%s:\n", identifier);

    printf("  n_sites  : %" PRId64 "\n", block.n_sites());
    if (block.sz_conserved() && block.charge_conserved()) {
      printf("  n_up     : %" PRId64 "\n", block.n_up());
      printf("  n_dn     : %" PRId64 "\n", block.n_dn());

    } else {
      printf("  n_up     : not conserved\n");
      printf("  n_dn     : not conserved\n");
    }

    if (block.symmetric()) {
      printf("  group    : defined with ID 0x%lx\n",
             (unsigned long)random::hash(block.permutation_group()));
      printf("  irrep    : defined with ID 0x%lx\n",
             (unsigned long)random::hash(block.irrep()));
    }

    std::stringstream ss;
    ss.imbue(std::locale("en_US.UTF-8"));
    ss << block.dim();
    printf("  dimension: %s\n", ss.str().c_str());

    std::stringstream ssmax;
    ssmax.imbue(std::locale("en_US.UTF-8"));
    ssmax << block.size_max();

    std::stringstream ssmin;
    ssmin.imbue(std::locale("en_US.UTF-8"));
    ssmin << block.size_min();

    std::stringstream ssavg;
    ssavg.imbue(std::locale("en_US.UTF-8"));
    ssavg << block.dim() / size;

    printf("  dimension (max local): %s\n", ssmax.str().c_str());
    printf("  dimension (min local): %s\n", ssmin.str().c_str());
    printf("  dimension (avg local): %s\n", ssavg.str().c_str());

    printf("  ID       : 0x%lx\n", (unsigned long)random::hash(block));
  }
}
#endif

void print_pretty(const char *identifier, Electron const &block) {
  printf("%s:\n", identifier);

  printf("  n_sites  : %" PRId64 "\n", block.n_sites());
  if (block.sz_conserved() && block.charge_conserved()) {
    printf("  n_up     : %" PRId64 "\n", block.n_up());
    printf("  n_dn     : %" PRId64 "\n", block.n_dn());

  } else {
    printf("  n_up     : not conserved\n");
    printf("  n_dn     : not conserved\n");
  }

  if (block.symmetric()) {
    printf("  group    : defined with ID 0x%lx\n",
           (unsigned long)random::hash(block.permutation_group()));
    printf("  irrep    : defined with ID 0x%lx\n",
           (unsigned long)random::hash(block.irrep()));
  }
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << block.size();
  printf("  dimension: %s\n", ss.str().c_str());
  printf("  ID       : 0x%lx\n", (unsigned long)random::hash(block));
}

void print_pretty(const char *identifier, RandomState const &rstate) {
  printf("%s:\n", identifier);
  printf("  RandomState, seed  : 0x%lx\n", (unsigned long)rstate.seed());
}

void print_pretty(const char *identifier, ProductState const &pstate) {
  printf("%s:\n", identifier);
  printf("  ProductState, n_sites: %" PRId64 "\n", pstate.n_sites());
  for (auto s : pstate) {
    printf("%s ", s.c_str());
  }
  printf("\n");
}

void print_pretty(const char *identifier, ProductState const &pstate);

void print_pretty(const char *identifier, Tmatrix const &tmat) {
  printf("%s:\n", identifier);
  print_pretty("alphas", tmat.alphas());
  print_pretty("betas", tmat.betas());
  print_pretty("eigenvalues", tmat.eigenvalues());
}

void print_pretty(const char *identifier, std::vector<int> const &vec) {
  printf("%s:\n", identifier);
  for (auto s : vec) {
    printf("%d ", s);
  }
  printf("\n");
}
void print_pretty(const char *identifier, std::vector<float> const &vec) {
  print_pretty(identifier, arma::Col<float>(vec));
}
void print_pretty(const char *identifier, std::vector<double> const &vec) {
  print_pretty(identifier, arma::Col<double>(vec));
}
void print_pretty(const char *identifier, std::vector<complex> const &vec) {
  print_pretty(identifier, arma::Col<complex>(vec));
}
void print_pretty(const char *identifier, std::vector<std::string> const &vec) {
  printf("%s:\n", identifier);
  for (auto s : vec) {
    printf("%s ", s.c_str());
  }
  printf("\n");
}

void print_pretty(const char *identifier, arma::uvec const &vec) {
  printf("%s:\n", identifier);
  vec.brief_print();
}

void print_pretty(const char *identifier, arma::ivec const &vec) {
  printf("%s:\n", identifier);
  vec.brief_print();
}

void print_pretty(const char *identifier, arma::umat const &mat) {
  printf("%s:\n", identifier);
  mat.brief_print();
}

void print_pretty(const char *identifier, arma::imat const &mat) {
  printf("%s:\n", identifier);
  mat.brief_print();
}
void print_pretty(const char *identifier, const arma::Mat<float> &mat) {
  printf("%s:\n", identifier);
  mat.brief_print();
}

void print_pretty(const char *identifier, const arma::Mat<double> &mat) {
  printf("%s:\n", identifier);
  mat.brief_print();
}

void print_pretty(const char *identifier,
                  const arma::Mat<std::complex<float>> &mat) {
  printf("%s:\n", identifier);
  mat.brief_print();
}

void print_pretty(const char *identifier,
                  const arma::Mat<std::complex<double>> &mat) {
  printf("%s:\n", identifier);
  mat.brief_print();
}

void print_pretty(const char *identifier, const arma::Col<float> &vec) {
  printf("%s:\n", identifier);
  vec.brief_print();
}

void print_pretty(const char *identifier, const arma::Col<double> &vec) {
  printf("%s:\n", identifier);
  vec.brief_print();
}

void print_pretty(const char *identifier,
                  const arma::Col<std::complex<float>> &vec) {
  printf("%s:\n", identifier);
  vec.brief_print();
}

void print_pretty(const char *identifier,
                  const arma::Col<std::complex<double>> &vec) {
  printf("%s:\n", identifier);
  vec.brief_print();
}

void print_pretty(const char *identifier, State const &state) {
  printf("%s:\n", identifier);
  if (state.isreal()) {
    printf("  real state\n");
  } else {
    printf("  cplx state\n");
  }
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << state.size();
  printf("  dimension: %s\n", ss.str().c_str());
  print_pretty((std::string(identifier) + std::string(".block()")).c_str(),
               state.block());
}

} // namespace xdiag::utils
