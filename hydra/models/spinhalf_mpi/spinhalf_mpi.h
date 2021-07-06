#pragma once

#include <vector>
#include <unordered_map>

#include <hydra/common.h>
#include <hydra/combinatorics/hashes.h>
#include <hydra/indexing/lintable.h>

namespace hydra {

template <class bit_t> class SpinhalfMPI {
public:
  SpinhalfMPI() = default;
  SpinhalfMPI(int n_sites, int n_up);

  inline int n_sites() const { return n_sites_; }

  inline bool sz_conserved() const { return sz_conserved_; }
  inline int sz() const { return sz_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }

  inline idx_t size() const { return size_; }
  inline idx_t dim() const { return dim_; }

  inline int process(bit_t prefix) { return hash_fnv1(prefix) % mpi_size_; }

  bool operator==(SpinhalfMPI const &rhs) const;
  bool operator!=(SpinhalfMPI const &rhs) const;

  int n_sites_;

  int n_prefix_bits_;
  int n_postfix_bits_;
  
  std::vector<bit_t> prefixes_;
  std::unordered_map<bit_t, std::pair<idx_t, idx_t>> prefix_limits_;
  std::vector<LinTable<bit_t>> postfix_lintables_;
  
private: 
  bool sz_conserved_;
  int n_up_;
  int n_dn_;
  int sz_;
 
  idx_t size_;
  idx_t dim_;
  
  int mpi_rank_;
  int mpi_size_;
};

} // namespace hydra
