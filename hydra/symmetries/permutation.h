#pragma once

#include <vector>
#include <hydra/io/toml/file_toml_handler.h>

namespace hydra {

class Permutation {
public:
  Permutation() = default;
  explicit Permutation(std::vector<int> const &array);
  explicit Permutation(std::initializer_list<int> list);
  explicit Permutation(io::FileTomlHandler && hdl);
  
  template <typename bit_t> bit_t apply(bit_t state) const;
  Permutation inverse() const;
  Permutation shuffle() const;

  int size() const;
  int operator[](int i) const;
  bool operator==(Permutation const &rhs) const;
  bool operator!=(Permutation const &rhs) const;
  
  std::vector<int> const& array() const;
private:
  int size_;
  std::vector<int> array_;
};

Permutation identity_permutation(int size);
Permutation operator*(Permutation const &p1, Permutation const &p2);
Permutation inverse(Permutation const &p);
Permutation shuffle(Permutation const &p);

} // namespace hydra
