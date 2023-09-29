#pragma once
#ifdef HYDRA_ENABLE_MPI

#include <hydra/common.h>

namespace hydra {

class tJDistributed {
public:
  using basis_t = basis_distributed_tj_variant_t;

  tJDistributed() = default;
  tJDistributed(int64_t n_sites, int64_t n_up, int64_t n_dn);

  int64_t n_sites() const;
  int64_t n_up() const;
  int64_t n_dn() const;

  int_t size() const;
  idx_t local_size() const;

  bool iscomplex(double precision = 1e-12) const;
  bool isreal(double precision = 1e-12) const;

  bool operator==(tJ const &rhs) const;
  bool operator!=(tJ const &rhs) const;
  basis_t const &basis() const;

private:
  int64_t n_sites_;
  int64_t n_up_;
  int64_t n_dn_;

  std::shared_ptr<basis_t> basis_;
  int64_t size_;
  int64_t local_size_;
}

} // namespace hydra

#endif
