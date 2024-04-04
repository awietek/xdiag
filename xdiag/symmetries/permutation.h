#pragma once

#include <vector>
#include <xdiag/io/toml/file_toml_handler.h>

namespace xdiag {

class Permutation {
public:
  Permutation() = default;
  explicit Permutation(std::vector<int64_t> const &array);
  explicit Permutation(std::initializer_list<int64_t> list);
  explicit Permutation(io::FileTomlHandler && hdl);
  
  template <typename bit_t> bit_t apply(bit_t state) const;
  Permutation inverse() const;
  Permutation shuffle() const;

  int64_t size() const;
  int64_t operator[](int64_t i) const;
  bool operator==(Permutation const &rhs) const;
  bool operator!=(Permutation const &rhs) const;
  
  std::vector<int64_t> const& array() const;
private:
  int64_t size_;
  std::vector<int64_t> array_;
};

Permutation identity_permutation(int64_t size);
Permutation operator*(Permutation const &p1, Permutation const &p2);
Permutation inverse(Permutation const &p);
Permutation shuffle(Permutation const &p);

} // namespace xdiag
