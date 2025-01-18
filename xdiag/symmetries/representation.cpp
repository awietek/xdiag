#include "representation.hpp"

#include <fstream>

#include <xdiag/utils/close.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

#include <iomanip>
#include <numeric>

namespace xdiag {

Representation::Representation(PermutationGroup const &group)
    : group_(group), characters_(arma::vec(group.size(), arma::fill::ones)) {}

template <typename T>
Representation::Representation(PermutationGroup const &group,
                               std::vector<T> const &characters)
    : group_(group), characters_(arma::Col<T>(characters)) {}

template Representation::Representation(PermutationGroup const &,
                                        std::vector<double> const &);
template Representation::Representation(PermutationGroup const &,
                                        std::vector<complex> const &);

template <typename T>
Representation::Representation(PermutationGroup const &group,
                               arma::Col<T> const &characters)
    : group_(group), characters_(characters) {}
template Representation::Representation(PermutationGroup const &,
                                        arma::vec const &);
template Representation::Representation(PermutationGroup const &,
                                        arma::cx_vec const &);

template <typename T>
Representation::Representation(PermutationGroup const &group, T *characters,
                               int64_t n_characters)
    : group_(group), characters_(arma::Col<T>(characters, n_characters, true)) {
}
template Representation::Representation(PermutationGroup const &group,
                                        double *characters,
                                        int64_t n_characters);
template Representation::Representation(PermutationGroup const &group,
                                        complex *characters,
                                        int64_t n_characters);

Representation::Representation(PermutationGroup const &group,
                               Vector const &characters)
    : group_(group), characters_(characters){};

Vector const &Representation::characters() const { return characters_; }

int64_t Representation::size() const { return characters_.size(); }
bool Representation::isreal() const { return characters_.is<arma::vec>(); }

bool Representation::operator==(Representation const &rhs) const {
  return (group_ == rhs.group_) && (characters_ == rhs.characters_);
}

bool Representation::operator!=(Representation const &rhs) const {
  return !operator==(rhs);
}

bool isreal(Representation const &irrep) { return irrep.isreal(); }

Representation trivial_representation(PermutationGroup const &group) {
  return Representation(group, arma::vec(group.size(), arma::fill::ones));
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
  return Representation();
}

Representation operator*(Representation const &r1,
                         Representation const &r2) try {
  return multiply(r1, r2);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return Representation();
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
