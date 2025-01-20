#include "basis_symmetric_np.hpp"

#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>

namespace xdiag::basis::tj {

template <typename bit_t>
template <typename coeff_t>
BasisSymmetricNp<bit_t>::BasisSymmetricNp(
    int64_t n_sites, int64_t nup, int64_t ndn, PermutationGroup const &group,
    arma::Col<coeff_t> const &characters) try
    : n_sites_(n_sites), n_up_(nup), n_dn_(ndn), group_action_(group),
      irrep_(group, characters),
      raw_ups_size_(combinatorics::binomial(n_sites, nup)),
      raw_dns_size_(combinatorics::binomial(n_sites, ndn)),
      raw_dnsc_size_(combinatorics::binomial(n_sites - nup, ndn)),
      lintable_ups_(n_sites, nup), lintable_dns_(n_sites, ndn),
      lintable_dnsc_(n_sites - nup, ndn), fermi_table_ups_(n_sites, nup, group),
      fermi_table_dns_(n_sites, ndn, group) {

  using combinatorics::Combinations;
  if ((nup + ndn) > n_sites) {
    XDIAG_THROW("nup + ndn > n_sites");
  } else if ((nup < 0) || (ndn < 0)) {
    XDIAG_THROW("nup < 0 or ndn < 0");
  } else if (n_sites < 0) {
    XDIAG_THROW("n_sites < 0");
  } else if (n_sites != group.n_sites()) {
    XDIAG_THROW("n_sites does not match the n_sites in PermutationGroup");
  }

  std::tie(reps_up_, idces_up_, syms_up_, sym_limits_up_) =
      symmetries::representatives_indices_symmetries_limits<bit_t>(
          lintable_ups_, group_action_);

  // if ups have trivial stabilizer, dns (compressed) are stored in front
  for (bit_t dns : Combinations<bit_t>(n_sites - nup, ndn)) {
    dns_storage_.push_back(dns);
    norms_storage_.push_back(1.0);
  }

  ups_offset_.resize(reps_up_.size());
  dns_limits_.resize(reps_up_.size());
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

  size_ = 0;
  int64_t idx_up = 0;

  for (bit_t ups : reps_up_) {
    ups_offset_[idx_up] = size_;
    auto syms = syms_ups(ups);

    // ups have trivial stabilizer -> dns stored in beginning
    if (syms.size() == 1) {
      span_size_t start = 0;
      span_size_t length = raw_dnsc_size_;
      dns_limits_[idx_up] = {start, length};
      size_ += length;
      // ups have non-trivial stabilizer, we store the dns configurations
    } else {
      bit_t not_ups = (~ups) & sitesmask;

      span_size_t start = dns_storage_.size();
      for (bit_t dnsc : Combinations<bit_t>(n_sites - nup, ndn)) {
        bit_t dns = bits::deposit(dnsc, not_ups);
        bit_t dns_rep =
            symmetries::representative_subset(dns, group_action_, syms);
        if (dns == dns_rep) {
          double norm = symmetries::norm_electron_subset(
              ups, dns, group_action_, characters, syms);
          if (norm > 1e-6) { // only keep dns with non-zero norm
            dns_storage_.push_back(dns_rep);
            norms_storage_.push_back(norm);
          }
        }
      }
      span_size_t end = dns_storage_.size();
      span_size_t length = end - start;
      dns_limits_[idx_up] = {start, length};
      size_ += length;
    }
    ++idx_up;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename bit_t>
inline int64_t BasisSymmetricNp<bit_t>::n_sites() const {
  return n_sites_;
}
template <typename bit_t> inline int64_t BasisSymmetricNp<bit_t>::n_up() const {
  return n_up_;
}
template <typename bit_t> inline int64_t BasisSymmetricNp<bit_t>::n_dn() const {
  return n_dn_;
}
template <typename bit_t>
GroupActionLookup<bit_t> const &BasisSymmetricNp<bit_t>::group_action() const {
  return group_action_;
}
template <typename bit_t>
Representation const &BasisSymmetricNp<bit_t>::irrep() const {
  return irrep_;
}
template <typename bit_t> int64_t BasisSymmetricNp<bit_t>::size() const {
  return size_;
}
template <typename bit_t> int64_t BasisSymmetricNp<bit_t>::dim() const {
  return size_;
}
template <typename bit_t>
typename BasisSymmetricNp<bit_t>::iterator_t
BasisSymmetricNp<bit_t>::begin() const {
  return iterator_t(*this, true);
}
template <typename bit_t>
typename BasisSymmetricNp<bit_t>::iterator_t
BasisSymmetricNp<bit_t>::end() const {
  return iterator_t(*this, false);
}
template <typename bit_t>
int64_t BasisSymmetricNp<bit_t>::index(bit_t ups, bit_t dns) const {
  auto syms = syms_ups(ups);
  int64_t idx_ups = index_ups(ups);
  int64_t up_offset = ups_offset(idx_ups);
  // trivial up-stabilizer (likely)
  if (syms.size() == 1) {
    bit_t sitesmask = ((bit_t)1 << n_sites_) - 1;
    bit_t not_ups = (~ups) & sitesmask;
    bit_t dnsc = bits::extract(dns, not_ups);
    int64_t idx_dns = dnsc_index(dnsc);
    return up_offset + idx_dns;
  }
  // non-trivial up-stabilizer (unlikely)
  else {
    auto dnss = dns_for_ups_rep(ups);
    auto [idx_dns, fermi_dn, sym] = index_dns_fermi_sym(dns, syms, dnss);
    return up_offset + idx_dns;
  }
}

template <typename bit_t> int64_t BasisSymmetricNp<bit_t>::n_rep_ups() const {
  return reps_up_.size();
}
template <typename bit_t>
bit_t BasisSymmetricNp<bit_t>::rep_ups(int64_t idx_up) const {
  return reps_up_[idx_up];
}
template <typename bit_t>
int64_t BasisSymmetricNp<bit_t>::ups_offset(int64_t idx_up) const {
  return ups_offset_[idx_up];
}

// index and fermi sign for dns with trivial stabilizer
template <typename bit_t>
std::pair<int64_t, bool>
BasisSymmetricNp<bit_t>::index_dns_fermi(bit_t dns, int64_t sym,
                                         bit_t not_ups) const {
  bit_t dns_rep = group_action_.apply(sym, dns);
  bit_t dns_rep_c = bits::extract(dns_rep, not_ups);
  int64_t idx_dns_rep = lintable_dnsc_.index(dns_rep_c);
  bool fermi_dns = fermi_table_dns_.sign(sym, dns);
  return {idx_dns_rep, fermi_dns};
}

template <typename bit_t>

std::pair<int64_t, bool>
BasisSymmetricNp<bit_t>::index_dns_fermi(bit_t dns, int64_t sym, bit_t not_ups,
                                         bit_t fermimask) const {
  bit_t dns_rep = group_action_.apply(sym, dns);
  bit_t dns_rep_c = bits::extract(dns_rep, not_ups);
  int64_t idx_dns_rep = lintable_dnsc_.index(dns_rep_c);
  bool fermi_dns = (bits::popcnt(dns & fermimask) & 1);
  fermi_dns ^= fermi_table_dns_.sign(sym, dns);
  return {idx_dns_rep, fermi_dns};
}

// index and fermi sign for dns with non-trivial stabilizer
template <typename bit_t>
std::tuple<int64_t, bool, int64_t> BasisSymmetricNp<bit_t>::index_dns_fermi_sym(
    bit_t dns, gsl::span<int64_t const> syms,
    gsl::span<bit_t const> dnss_out) const {
  auto [rep_dns, rep_sym] =
      symmetries::representative_sym_subset(dns, group_action_, syms);
  auto it = std::lower_bound(dnss_out.begin(), dnss_out.end(), rep_dns);
  if ((it != dnss_out.end()) && (*it == rep_dns)) {
    bool fermi_dns = fermi_table_dns_.sign(rep_sym, dns);
    return {std::distance(dnss_out.begin(), it), fermi_dns, rep_sym};
  } else {
    return {invalid_index, false, rep_sym};
  }
}

template <typename bit_t>

std::tuple<int64_t, bool, int64_t> BasisSymmetricNp<bit_t>::index_dns_fermi_sym(
    bit_t dns, gsl::span<int64_t const> syms, gsl::span<bit_t const> dnss_out,
    bit_t fermimask) const {
  auto [rep_dns, rep_sym] =
      symmetries::representative_sym_subset(dns, group_action_, syms);
  auto it = std::lower_bound(dnss_out.begin(), dnss_out.end(), rep_dns);
  if ((it != dnss_out.end()) && (*it == rep_dns)) {
    bool fermi_dns = (bits::popcnt(dns & fermimask) & 1);
    fermi_dns ^= fermi_table_dns_.sign(rep_sym, dns);
    return {std::distance(dnss_out.begin(), it), fermi_dns, rep_sym};
  } else {
    return {invalid_index, false, rep_sym};
  }
}

// Retrieving index of representative and symmetries for ups
template <typename bit_t>

int64_t BasisSymmetricNp<bit_t>::index_ups(bit_t ups) const {
  return idces_up_[lintable_ups_.index(ups)];
}

template <typename bit_t>
gsl::span<int64_t const> BasisSymmetricNp<bit_t>::syms_ups(bit_t ups) const {
  int64_t idx_ups = lintable_ups_.index(ups);
  auto [start, length] = sym_limits_up_[idx_ups];
  return {syms_up_.data() + start, length};
}

template <typename bit_t>
std::pair<int64_t, gsl::span<int64_t const>>
BasisSymmetricNp<bit_t>::index_syms_up(bit_t ups) const {
  int64_t idx_ups = lintable_ups_.index(ups);
  int64_t index = idces_up_[idx_ups];
  auto [start, length] = sym_limits_up_[idx_ups];
  return {index, {syms_up_.data() + start, length}};
}

// Retrieving dns states and norms for given up configuration
template <typename bit_t>
gsl::span<bit_t const>
BasisSymmetricNp<bit_t>::dns_for_ups_rep(bit_t ups) const {
  int64_t idx_ups = index_ups(ups);
  auto [start, length] = dns_limits_[idx_ups];
  return {dns_storage_.data() + start, length};
}

template <typename bit_t>
gsl::span<double const>
BasisSymmetricNp<bit_t>::norms_for_ups_rep(bit_t ups) const {
  int64_t idx_ups = index_ups(ups);
  auto [start, length] = dns_limits_[idx_ups];
  return {norms_storage_.data() + start, length};
}

template <typename bit_t>
std::pair<gsl::span<bit_t const>, gsl::span<double const>>
BasisSymmetricNp<bit_t>::dns_norms_for_up_rep(bit_t ups) const {
  int64_t idx_ups = index_ups(ups);
  auto [start, length] = dns_limits_[idx_ups];
  auto dnss = gsl::span<bit_t const>{dns_storage_.data() + start, length};
  auto norms = gsl::span<double const>{norms_storage_.data() + start, length};
  return {dnss, norms};
}

// Fermi sign when applying sym on states
template <typename bit_t>

bool BasisSymmetricNp<bit_t>::fermi_bool_ups(int64_t sym, bit_t ups) const {
  return fermi_table_ups_.sign(sym, ups);
}
template <typename bit_t>

bool BasisSymmetricNp<bit_t>::fermi_bool_dns(int64_t sym, bit_t dns) const {
  return fermi_table_dns_.sign(sym, dns);
}

template <typename bit_t>
int64_t BasisSymmetricNp<bit_t>::dnsc_index(bit_t dns) const {
  return lintable_dnsc_.index(dns);
}

template class BasisSymmetricNp<uint32_t>;
template class BasisSymmetricNp<uint64_t>;
template BasisSymmetricNp<uint32_t>::BasisSymmetricNp<double>(
    int64_t, int64_t, int64_t, PermutationGroup const &,
    arma::Col<double> const &);
template BasisSymmetricNp<uint32_t>::BasisSymmetricNp<complex>(
    int64_t, int64_t, int64_t, PermutationGroup const &,
    arma::Col<complex> const &);
template BasisSymmetricNp<uint64_t>::BasisSymmetricNp<double>(
    int64_t, int64_t, int64_t, PermutationGroup const &,
    arma::Col<double> const &);
template BasisSymmetricNp<uint64_t>::BasisSymmetricNp<complex>(
    int64_t, int64_t, int64_t, PermutationGroup const &,
    arma::Col<complex> const &);

template <typename bit_t>
BasisSymmetricNpIterator<bit_t>::BasisSymmetricNpIterator(
    BasisSymmetricNp<bit_t> const &basis, bool begin)
    : basis_(basis), sitesmask_(((bit_t)1 << basis.n_sites()) - 1),
      up_idx_(begin ? 0 : basis.n_rep_ups()), dn_idx_(0) {
  if (basis.dim() == 0) {
    up_idx_ = 0;
  }

  if ((basis.n_rep_ups() > 0) && begin) {
    dns_for_ups_rep_ = basis.dns_for_ups_rep(basis.rep_ups(0));
  }
}

template <typename bit_t>
BasisSymmetricNpIterator<bit_t> &BasisSymmetricNpIterator<bit_t>::operator++() {
  ++dn_idx_;
  if (dn_idx_ == dns_for_ups_rep_.size()) {
    dn_idx_ = 0;
    do {
      ++up_idx_;
      if (up_idx_ == basis_.n_rep_ups()) {
        return *this;
      }
      bit_t ups = basis_.rep_ups(up_idx_);
      dns_for_ups_rep_ = basis_.dns_for_ups_rep(ups);
    } while (dns_for_ups_rep_.size() == 0);
  }
  return *this;
}

template <typename bit_t>
std::pair<bit_t, bit_t> BasisSymmetricNpIterator<bit_t>::operator*() const {
  bit_t ups = basis_.rep_ups(up_idx_);
  auto syms = basis_.syms_ups(ups);
  if (syms.size() == 1) {
    bit_t not_ups = (~ups) & sitesmask_;
    bit_t dnsc = dns_for_ups_rep_[dn_idx_];
    bit_t dns = bits::deposit(dnsc, not_ups);
    return {ups, dns};
  } else {
    bit_t dns = dns_for_ups_rep_[dn_idx_];
    return {ups, dns};
  }
}

template <typename bit_t>
bool BasisSymmetricNpIterator<bit_t>::operator!=(
    BasisSymmetricNpIterator<bit_t> const &rhs) const {
  return (up_idx_ != rhs.up_idx_) || (dn_idx_ != rhs.dn_idx_);
}

template class BasisSymmetricNpIterator<uint32_t>;
template class BasisSymmetricNpIterator<uint64_t>;

} // namespace xdiag::basis::tj
