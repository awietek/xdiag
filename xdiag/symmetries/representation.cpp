#include "representation.hpp"

#include <fstream>

#include <xdiag/utils/close.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

#include <numeric>
#include <iomanip>

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

bool Representation::operator==(Representation const &rhs) const {
  return close(arma::cx_vec(characters_), arma::cx_vec(rhs.characters_));
}

bool Representation::operator!=(Representation const &rhs) const {
  return !operator==(rhs);
}

Representation read_representation(std::string filename, std::string repname) {
  std::ifstream File(filename.c_str());
  if (File.fail()) {
    std::cerr << "Error in read_charactertable: Could not open "
              << "file with filename [" << filename << "] given. Abort."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  std::vector<std::string> names;
  std::vector<std::vector<int64_t>> allowed_symmetries_arr;
  std::vector<std::vector<complex>> characters_arr;

  std::string tobeparsed;
  std::string::size_type pos;

  // Jump to Irreps and parse nreps
  File >> tobeparsed;
  while (tobeparsed.find("[Irreps]") == std::string::npos)
    File >> tobeparsed;
  pos = tobeparsed.find('=');
  int64_t nreps;
  if (pos != std::string::npos)
    nreps = atoi(tobeparsed.substr(pos + 1, std::string::npos).c_str());
  else
    nreps = -1;
  assert(nreps >= 0);

  // Loop over all representations
  bool found = false;
  for (int64_t i = 0; i < nreps; ++i) {
    File >> tobeparsed;
    while (tobeparsed.find("[Representation]") == std::string::npos)
      File >> tobeparsed;
    pos = tobeparsed.find('=');

    // Get name of representation
    std::string name = tobeparsed.substr(pos + 1, std::string::npos);
    std::vector<int64_t> allowed_ops;
    std::vector<complex> characters;

    // parse number of allowed operations
    while (tobeparsed.find("[AllowedOps]") == std::string::npos)
      File >> tobeparsed;
    pos = tobeparsed.find('=');
    int64_t n_allowed_ops;
    if (pos != std::string::npos)
      n_allowed_ops =
          atoi(tobeparsed.substr(pos + 1, std::string::npos).c_str());
    else
      n_allowed_ops = -1;

    // parse allowed operations
    for (int64_t so = 0; so < n_allowed_ops; ++so) {
      int64_t w;
      File >> w;
      if (!File.good()) {
        std::cerr << "Read Error in read_representation (I)" << std::endl;
        exit(EXIT_FAILURE);
      }
      allowed_ops.push_back(w);
    }

    if (!File.good()) {
      std::cerr << "Read Error in read_representation (II)" << std::endl;
      exit(EXIT_FAILURE);
    }

    // parse bloch factors
    for (int64_t so = 0; so < n_allowed_ops; ++so) {
      double re, im;
      File >> re >> im;
      if (!File.good()) {
        std::cerr << "Read Error in read_representation (III)" << std::endl;
        exit(EXIT_FAILURE);
      }
      characters.push_back(re + complex(0, 1) * im);
    }

    auto rep = Representation(characters, allowed_ops);
    if (name == repname) {
      found = true;
      return rep;
    }
  }
  if (!found)
    Log.err("Error reading representations: "
            "name not found in file");

  return Representation();
}

bool Representation::isreal(double precision) const {
  for (auto c : characters_) {
    if (std::abs(imag(c)) > precision) {
      return false;
    }
  }
  return true;
}

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
