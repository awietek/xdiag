#pragma once

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/linalg/lanczos/tmatrix.h>
#include <lila/all.h>

namespace hydra::utils {

  void PrintPretty(const char* identifier, Bond const& bond);
  void PrintPretty(const char* identifier, BondList const& bondlist);
  void PrintPretty(const char* identifier, Couplings const& couplings);

  template <typename T>
  void PrintPretty(const char* identifier, lila::Vector<T> const& vector){
    lila::PrintPretty(identifier, vector);
  }
  template <typename T>
  void PrintPretty(const char* identifier, lila::Matrix<T> const& matrix){
    lila::PrintPretty(identifier, matrix);
  }

  void PrintPretty(const char* identifier, Tmatrix const& tmat);
  
} // namespace hydra::utils
