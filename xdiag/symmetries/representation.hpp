#pragma once

#include <string>
#include <vector>

#include <xdiag/common.hpp>
#include <xdiag/io/toml/file_toml_handler.hpp>
#include <xdiag/symmetries/permutation_group.hpp>

namespace xdiag {

class Representation {
public:
  Representation() = default;
  Representation(std::initializer_list<complex> list);
  explicit Representation(std::vector<complex> const &characters);
  Representation(std::vector<complex> const &characters,
                 std::vector<int64_t> const &allowed_symmetries);
  explicit Representation(io::FileTomlHandler &&hdl);
  Representation(complex const *characters, int64_t n_characters); 

  complex character(int64_t idx) const;
  std::vector<int64_t> allowed_symmetries() const;
  std::vector<complex> const &characters() const;
  std::vector<double> const &characters_real() const;

  Representation subgroup(std::vector<int64_t> const &symmetry_numbers) const;
  int64_t size() const;

  bool isreal(double precision = 1e-12) const;

  bool operator==(Representation const &rhs) const;
  bool operator!=(Representation const &rhs) const;

private:
  std::vector<complex> characters_;
  std::vector<double> characters_real_;
  std::vector<int64_t> allowed_symmetries_;
};

Representation trivial_representation(int64_t size);
Representation trivial_representation(PermutationGroup const &group);

PermutationGroup allowed_subgroup(PermutationGroup const &group,
                                  Representation const &irrep);

Representation multiply(Representation const &r1, Representation const &r2);
Representation operator*(Representation const &r1, Representation const &r2);

std::ostream &operator<<(std::ostream &out, Representation const &irrep);
std::string to_string(Representation const &irrep);
} // namespace xdiag
