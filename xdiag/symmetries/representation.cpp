#include "representation.hpp"

#include <fstream>

#include <xdiag/utils/close.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

#include <iomanip>
#include <numeric>

namespace xdiag {

Representation::Representation(std::initializer_list<complex> list)
    : Representation(std::vector<complex>(list)) {}

Representation::Representation(std::vector<complex> const &characters)
    : characters_(characters), characters_real_(characters.size()),
      allowed_symmetries_(characters.size(), 0) {
  for (int64_t idx = 0; idx < (int64_t)characters.size(); ++idx) {
    characters_real_[idx] = std::real(characters[idx]);
  }
  std::iota(allowed_symmetries_.begin(), allowed_symmetries_.end(), 0);
}

Representation::Representation(complex const *characters, int64_t n_characters)
    : Representation(
          std::vector<complex>(characters, characters + n_characters)) {}

Representation::Representation(std::vector<complex> const &characters,
                               std::vector<int64_t> const &allowed_symmetries)
    : characters_(characters), characters_real_(characters.size()),
      allowed_symmetries_(allowed_symmetries) {

  assert(characters.size() == allowed_symmetries.size());

  for (int64_t idx = 0; idx < (int64_t)characters.size(); ++idx) {
    characters_real_[idx] = std::real(characters[idx]);
  }
}

Representation::Representation(io::FileTomlHandler &&hdl)
    : Representation(hdl.as<Representation>()) {}

complex Representation::character(int64_t idx) const {
  return characters_.at(idx);
}
std::vector<int64_t> Representation::allowed_symmetries() const {
  return allowed_symmetries_;
}

std::vector<complex> const &Representation::characters() const {
  return characters_;
}
std::vector<double> const &Representation::characters_real() const {
  return characters_real_;
}

Representation
Representation::subgroup(std::vector<int64_t> const &symmetry_numbers) const {
  std::vector<complex> sub_characters;
  for (auto sym : symmetry_numbers)
    sub_characters.push_back(characters_[sym]);
  return Representation(sub_characters);
}
int64_t Representation::size() const { return characters_.size(); }
bool Representation::isreal() const { return characters_.is<arma::vec>(); }

bool Representation::operator==(Representation const &rhs) const {
  return close(arma::cx_vec(characters_), arma::cx_vec(rhs.characters_));
}

bool Representation::operator!=(Representation const &rhs) const {
  return !operator==(rhs);
}

bool isreal(Representation const &irrep) { return irrep.isreal(); }

Representation trivial_representation(int64_t size) {
  return Representation(std::vector<complex>(size, {1.0, 0.0}));
}

Representation trivial_representation(PermutationGroup const &group) {
  return Representation(std::vector<complex>(group.size(), {1.0, 0.0}));
}

Representation multiply(Representation const &r1,
                        Representation const &r2) try {
  if (r1.size() != r2.size()) {
    XDIAG_THROW("Cannot construct product Representation, sizes of two inputs "
                "are not equal.");
  }

  auto c1 = r1.characters();
  auto c2 = r2.characters();
  auto c3 = std::vector<complex>(r1.size());

  for (int64_t i = 0; i < (int64_t)r1.size(); ++i) {
    c3[i] = c1[i] * c2[i];
  }

  return Representation(c3);
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

PermutationGroup allowed_subgroup(PermutationGroup const &group,
                                  Representation const &irrep) {

  auto const &allowed_symmetries = irrep.allowed_symmetries();

  if (allowed_symmetries.size() > 0) {
    return group.subgroup(allowed_symmetries);
  } else {
    return group;
  }
}

std::ostream &operator<<(std::ostream &out, Representation const &irrep) {
  out << "size      : " << irrep.size() << "\n";
  out << "characters:\n";
  for (auto c : irrep.characters()) {
    out << std::setprecision(8) << c << "\n";
  }
  return out;
}
std::string to_string(Representation const &irrep) {
  return to_string_generic(irrep);
}

} // namespace xdiag
