//=======================================================================
//                 Alea library $Revision: 1.1.1.1 $
//          (c) 1994-1998 by Matthias Troyer and Beat Ammon
//                ISSP, University of Tokyo, Japan
//		  CSCS, ETH Zurich, Switzerland
//
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Institutions (University of Tokyo, ETH Zurich)
// nor the Authors make any representations about the suitability of this
// software for any purpose.  This software is provided ``as is'' without
// express or implied warranty.
// 
//=======================================================================

#ifndef HYDRA_PARAMETERS_OSIRIS_PARSE_
#define HYDRA_PARAMETERS_OSIRIS_PARSE_

#include <iostream>
#include <complex>
#include "parameters.h"

//=======================================================================
// parser object
//
// implements a parser on an istream
//-----------------------------------------------------------------------

namespace hydra { namespace parameters {

class parser {
public:
  using token = char;
private:
  std::istream& in;
  parameter_value val;
  char c;
  
  char put_back_token;
  bool was_put_back;

  void parse_string();
  char parse_ident();
  char parse_number();
  char parse_complex();
  void parse_collection();

  bool isskip(char);
  char next_char();
  char eat_ws();

public:
  parser(std::istream& i) : in(i), was_put_back(false) { in >> c;}
  char next_token();
  char next_token_nows();
  bool putback(char);
  const parameter_value& value() const {return val;}
};


}} // namespace hydra::parameters

#endif
