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

#pragma once

//=======================================================================
// This file includes implementation definitions for the parameter
// classes
//=======================================================================

#include "parameters.h"
#include "parser.h"

//=======================================================================
// parameters_output
//
// is a function object to output parameters
//-----------------------------------------------------------------------
namespace hydra {

class parameters_output {
private:
  std::ostream &out;

public:
  parameters_output(std::ostream &o) : out(o) {}
  inline void
  operator()(const std::pair<const std::string, parameter_value> &x) const {
    out << x.first << " = " << x.second << ";\n";
  }
};

class parameters_collection_output {
private:
  std::ostream &out;

public:
  parameters_collection_output(std::ostream &o) : out(o) {}
  inline void
  operator()(const std::pair<const std::string, Parameters> &x) const {
    out << x.first << " {\n";
    out << x.second;
    out << "} " << x.first << "\n";
  }
};

//=======================================================================
// functions for reading parameters from a parser
//-----------------------------------------------------------------------

parser &operator>>(parser &, Parameters &);
parser &operator>>(parser &, parameters_collection &);

} // namespace hydra
