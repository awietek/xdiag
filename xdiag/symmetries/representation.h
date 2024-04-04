#pragma once

#include <string>
#include <vector>

#include <xdiag/common.h>
#include <xdiag/io/toml/file_toml_handler.h>
#include <xdiag/symmetries/permutation_group.h>

namespace xdiag {

class Representation {
public:
  Representation() = default;
  explicit Representation(std::vector<complex> const &characters);
  Representation(std::vector<complex> const &characters,
                 std::vector<int64_t> const &allowed_symmetries);
  explicit Representation(io::FileTomlHandler &&hdl);

  complex character(int64_t idx) const { return characters_.at(idx); }
  std::vector<int64_t> allowed_symmetries() const { return allowed_symmetries_; }

  std::vector<complex> const &characters() const { return characters_; }
  std::vector<double> const &characters_real() const {
    return characters_real_;
  }
  Representation subgroup(std::vector<int64_t> const &symmetry_numbers) const;

  int64_t size() const { return characters_.size(); }

  bool iscomplex(double precision = 1e-12) const;
  bool isreal(double precision = 1e-12) const;

  bool operator==(Representation const &rhs) const;
  bool operator!=(Representation const &rhs) const;

private:
  std::vector<complex> characters_;
  std::vector<double> characters_real_;
  std::vector<int64_t> allowed_symmetries_;
};

Representation read_representation(std::string filename, std::string repname);
Representation trivial_representation(int64_t size);
Representation trivial_representation(PermutationGroup const &group);

PermutationGroup allowed_subgroup(PermutationGroup const &group,
                                  Representation const &irrep);

Representation operator*(Representation const &r1, Representation const &r2);

} // namespace xdiag
