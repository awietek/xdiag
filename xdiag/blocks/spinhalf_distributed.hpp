#pragma once
#ifdef XDIAG_USE_MPI
#include <xdiag/basis/spinhalf_distributed/basis_spinhalf_distributed.hpp>
#include <xdiag/common.hpp>

namespace xdiag {

class SpinhalfDistributed {
public:
  using basis_t = basis::BasisSpinhalfDistributed;
  SpinhalfDistributed() = default;
  SpinhalfDistributed(int64_t n_sites, int64_t nup);

  int64_t n_sites() const;
  int64_t n_up() const;

  int64_t dim() const;
  int64_t size() const;
  int64_t size_max() const;
  int64_t size_min() const;

  bool isreal(double precision = 1e-12) const;

  bool operator==(SpinhalfDistributed const &rhs) const;
  bool operator!=(SpinhalfDistributed const &rhs) const;

  basis_t const &basis() const;

private:
  int64_t n_sites_;
  int64_t n_up_;

  std::shared_ptr<basis_t> basis_;
  int64_t dim_;
  int64_t size_;
};

} // namespace xdiag
#endif
