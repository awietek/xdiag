#include "representation.hpp"

#include <fstream>
#include <iomanip>
#include <numeric>

namespace xdiag {

template <typename T>
static void check_characters(PermutationGroup const &group,
                             arma::Col<T> const &characters,
                             double precision = 1e-12) try {
  int64_t n = group.size();
  for (int64_t i = 0; i < n; ++i) {
    for (int64_t j = 0; j < n; ++j) {
      T cij = characters[group.multiply(i, j)];
      T cicj = characters[i] * characters[j];
      if (std::abs(cij - cicj) > precision) {
        XDIAG_THROW(
            "Characters of Representation are not fulfilling multiplicative "
            "law with respect to the PermutationGroup, i.e. there exists group "
            "elements \"g\" and \"h\" such that c(g)*c(h) != c(gh)")
      }
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Representation::Representation(PermutationGroup const &group,
                               Vector const &characters) try
    : group_(group), characters_(characters) {
  if (characters.isreal()) {
    check_characters(group, characters.as<arma::vec>());
  } else {
    check_characters(group, characters.as<arma::cx_vec>());

    // If imaginary part is small, make it real
    auto ni = arma::norm(arma::imag(characters.as<arma::cx_vec>()));
    if (ni < 1e-14) {
      characters_ = Vector(arma::real(characters.as<arma::cx_vec>()));
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Representation::Representation(PermutationGroup const &group) try
    : Representation(group, Vector(arma::vec(group.size(), arma::fill::ones))) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename T>
Representation::Representation(PermutationGroup const &group,
                               std::vector<T> const &characters) try
    : Representation(group, Vector(arma::Col<T>(characters))) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template Representation::Representation(PermutationGroup const &,
                                        std::vector<double> const &);
template Representation::Representation(PermutationGroup const &,
                                        std::vector<complex> const &);

template <typename T>
Representation::Representation(PermutationGroup const &group,
                               arma::Col<T> const &characters) try
    : Representation(group, Vector(characters)) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template Representation::Representation(PermutationGroup const &,
                                        arma::vec const &);
template Representation::Representation(PermutationGroup const &,
                                        arma::cx_vec const &);

template <typename T>
Representation::Representation(PermutationGroup const &group, T *characters,
                               int64_t n_characters) try
    : Representation(group,
                     Vector(arma::Col<T>(characters, n_characters, true))) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template Representation::Representation(PermutationGroup const &group,
                                        double *characters,
                                        int64_t n_characters);
template Representation::Representation(PermutationGroup const &group,
                                        complex *characters,
                                        int64_t n_characters);
int64_t Representation::size() const { return characters_.size(); }
bool Representation::operator==(Representation const &rhs) const {
  return (group_ == rhs.group_) && (characters_ == rhs.characters_);
}

bool Representation::operator!=(Representation const &rhs) const {
  return !operator==(rhs);
}

PermutationGroup const &Representation::group() const { return group_; }
Vector const &Representation::characters() const { return characters_; }
bool Representation::isreal() const { return characters_.is<arma::vec>(); }

int64_t size(Representation const &irrep) { return irrep.size(); }
bool isreal(Representation const &irrep) { return irrep.isreal(); }
bool isapprox(Representation const &r1, Representation const &r2, double rtol,
              double atol) {
  return (r1.group() == r2.group()) &&
         isapprox(r1.characters(), r2.characters(), rtol, atol);
}
Representation multiply(Representation const &r1,
                        Representation const &r2) try {
  if (r1.group() != r2.group()) {
    XDIAG_THROW(
        "The two given representations are not defined for the same group.");
  }

  if (isreal(r1) && isreal(r2)) {
    auto c1 = r1.characters().as<arma::vec>();
    auto c2 = r2.characters().as<arma::vec>();
    return Representation(
        r1.group(),
        arma::vec(c1 % c2)); // % means element wise multiplication in armadillo
  } else {
    auto c1 = r1.characters().as<arma::cx_vec>();
    auto c2 = r2.characters().as<arma::cx_vec>();
    return Representation(r1.group(), arma::cx_vec(c1 % c2));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Representation operator*(Representation const &r1,
                         Representation const &r2) try {
  return multiply(r1, r2);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::ostream &operator<<(std::ostream &out, Representation const &irrep) {
  out << "size      : " << irrep.size() << "\n";
  out << "characters:\n";

  if (isreal(irrep)) {
    for (auto c : irrep.characters().as<arma::vec>()) {
      out << std::setprecision(8) << c << "\n";
    }
  } else {
    for (auto c : irrep.characters().as<arma::cx_vec>()) {
      out << std::setprecision(8) << c << "\n";
    }
  }
  return out;
}
std::string to_string(Representation const &irrep) {
  return to_string_generic(irrep);
}

} // namespace xdiag
