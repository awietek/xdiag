#pragma once

#include <ostream>
#include <string>
#include <vector>

#include <hydra/operators/bond_types.h>
#include <hydra/parameters/parameters.h>

namespace hydra {

class Bond {
public:
  Bond(std::string type, std::string coupling, std::vector<int> const &sites);
  Bond(std::string type, std::string coupling, std::vector<int> const &sites,
       Parameters const &parameters);
  inline std::string type() const { return type_; }
  inline std::string coupling() const { return coupling_; }
  inline std::vector<int> sites() const { return sites_; }
  inline int site(int j) const { return sites_.at(j); }
  inline int size() const { return (int)sites_.size(); }
  inline bool has_parameters() const { return has_parameters_; }
  inline Parameters parameters() const { return parameters_; }
  inline int operator[](int j) const { return site(j); }

private:
  std::string type_;
  std::string coupling_;
  std::vector<int> sites_;
  bool has_parameters_;
  Parameters parameters_;
};

std::vector<int> common_sites(Bond b1, Bond b2);

std::ostream &operator<<(std::ostream &out, const Bond &bond);

bool operator==(const Bond &lhs, const Bond &rhs);

struct TypeCoupling {
  TypeCoupling(const std::string &type, const std::string &coupling);
  inline std::string type() const { return type_; };
  inline std::string coupling() const { return coupling_; };

  std::string type_;
  std::string coupling_;
};

std::ostream &operator<<(std::ostream &out, const TypeCoupling &tc);

bool operator==(const TypeCoupling &tc1, const TypeCoupling &tc2);

TypeCoupling type_coupling(Bond const &bond);

bool is_complex(Bond const &bond);
bool is_real(Bond const &bond);

} // namespace hydra