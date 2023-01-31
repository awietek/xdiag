#pragma once

#include <string>
#include <vector>

#include <hydra/common.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/io/toml/file_toml_handler.h>

namespace hydra {

class Representation {
public:
  Representation() = default;
  explicit Representation(std::vector<complex> const &characters);
  Representation(std::vector<complex> const &characters,
                 std::vector<int> const &allowed_symmetries);
  explicit Representation(io::FileTomlHandler && hdl);

  complex character(int idx) const { return characters_.at(idx); }
  std::vector<int> allowed_symmetries() const { return allowed_symmetries_; }

  std::vector<complex> const &characters() const { return characters_; }
  std::vector<double> const &characters_real() const {
    return characters_real_;
  }
  Representation subgroup(std::vector<int> const &symmetry_numbers) const;

  idx_t size() const { return characters_.size(); }

  bool operator==(Representation const &rhs) const;
  bool operator!=(Representation const &rhs) const;

private:
  std::vector<complex> characters_;
  std::vector<double> characters_real_;
  std::vector<int> allowed_symmetries_;
};

Representation read_representation(std::string filename, std::string repname);

bool is_complex(Representation const &cpls);
bool is_real(Representation const &cpls);

Representation trivial_representation(idx_t size);
Representation trivial_representation(PermutationGroup const& group);

PermutationGroup allowed_subgroup(PermutationGroup const &group,
                                  Representation const &irrep);
  
Representation operator*(Representation const &r1, Representation const &r2);

} // namespace hydra
