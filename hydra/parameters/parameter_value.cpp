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

//=======================================================================
// This file includes the class definitions for the parameter classes
//=======================================================================

#include "parameters.h"

#include <iomanip>
#include <hydra/utils/logger.h>

//=======================================================================
// parameter_value
//
// a parameter value, containing the variable type and value
//-----------------------------------------------------------------------

namespace hydra {

//-----------------------------------------------------------------------
// CONSTRUCTORS AND DESTRUCTOR
//-----------------------------------------------------------------------

// invalid parameter

parameter_value::parameter_value() : the_type(p_is_invalid) {}

// make a string parameter

parameter_value::parameter_value(const std::string &s)
    : the_string(s), the_type(p_is_string) {}

// make an integer parameter

parameter_value::parameter_value(const long l)
    : valint(l), the_type(p_is_integer) {}

// make a floating point parameter

parameter_value::parameter_value(const double d)
    : valfloat(d), the_type(p_is_float) {}

// make a complex parameter

parameter_value::parameter_value(const complex &d)
    : valcomplex(d), the_type(p_is_complex) {}

// make a bool parameter

parameter_value::parameter_value(const bool b)
    : valbool(b), the_type(p_is_bool) {}

// copy a parameter

parameter_value::parameter_value(const parameter_value &p)
    : the_string(p.the_string), the_type(p.the_type) {
  switch (the_type) {
  case p_is_integer:
    valint = p.valint;
    break;

  case p_is_float:
    valfloat = p.valfloat;
    break;

  case p_is_string:
    break;

  case p_is_invalid:
    break;

  case p_is_complex:
    valcomplex = p.valcomplex;
    break;

  case p_is_bool:
    valbool = p.valbool;
    break;

  default:
    Log.err("illegal parameter type in parameter_value::parameter_value");
  }
}

//-----------------------------------------------------------------------
// ASSIGNMENT OPERATORS
//-----------------------------------------------------------------------

// assign another parameter

parameter_value &parameter_value::operator=(const parameter_value &p) {
  the_string = p.the_string;
  the_type = p.the_type;

  switch (the_type) {
  case p_is_integer:
    valint = p.valint;
    break;

  case p_is_float:
    valfloat = p.valfloat;
    break;

  case p_is_string:
    break;

  case p_is_invalid:
    break;

  case p_is_complex:
    valcomplex = p.valcomplex;
    break;

  case p_is_bool:
    valbool = p.valbool;
    break;

  default:
    Log.err("illegal parameter type in parameter_value::operator=");
  }

  return *this;
}

// assign a string

parameter_value &parameter_value::operator=(const std::string &s) {
  the_type = p_is_string;
  the_string = s;
  return *this;
}

parameter_value &parameter_value::operator=(const char *s) {
  the_type = p_is_string;
  the_string = s;
  return *this;
}

// assign a floating point number

parameter_value &parameter_value::operator=(const double s) {
  the_type = p_is_float;
  valfloat = s;
  return *this;
}

// assign an integer

parameter_value &parameter_value::operator=(const long s) {
  the_type = p_is_integer;
  valint = s;
  return *this;
}

// assign a bool

parameter_value &parameter_value::operator=(const bool b) {
  the_type = p_is_bool;
  valbool = b;
  return *this;
}

// assign a complex number

parameter_value &parameter_value::operator=(const complex &c) {
  the_type = p_is_complex;
  valcomplex = c;
  return *this;
}

//-----------------------------------------------------------------------
// CONVERSION OPERATORS
//-----------------------------------------------------------------------

// convert to a character array

parameter_value::operator const char *() const {
  if (the_type != p_is_string)
    type_error(p_is_string);

  return the_string.c_str();
}

// convert to a string

parameter_value::operator const std::string &() const {
  if (the_type != p_is_string)
    type_error(p_is_string);

  return the_string;
}

// convert to a long integer

parameter_value::operator long() const {
  if (the_type == p_is_integer)
    return valint;
  else if (the_type == p_is_bool)
    return valbool;
  else
    type_error(p_is_integer);
  return 0; // dummy;
}

// convert to integer

parameter_value::operator int() const {
  if (the_type == p_is_integer)
    return valint;
  else if (the_type == p_is_bool)
    return valbool;
  else
    type_error(p_is_integer);
  return 0; // dummy
}

// convert to floating point number

parameter_value::operator float() const {
  if (the_type == p_is_float)
    return valfloat;
  else if (the_type == p_is_integer)
    return valint;
  else if (the_type == p_is_bool)
    return valbool;
  else
    type_error(p_is_float);
  return 0; // dummy
}

// convert to double precision floating point number

parameter_value::operator double() const {
  if (the_type == p_is_float)
    return valfloat;
  else if (the_type == p_is_integer)
    return valint;
  else if (the_type == p_is_bool)
    return valbool;
  else
    type_error(p_is_float);
  return 0; // dummy
}

// convert to complex

parameter_value::operator complex() const {
  if (the_type == p_is_complex)
    return valcomplex;
  else if (the_type == p_is_float)
    return complex(valfloat, 0.);
  else if (the_type == p_is_integer)
    return complex(valint, 0.);
  else if (the_type == p_is_bool)
    return complex(double(valbool), 0.);
  else
    type_error(p_is_complex);
  return complex(0, 0); // dummy
}

// convert to bool

parameter_value::operator bool() const {
  if (the_type == p_is_integer)
    return valint;
  else if (the_type == p_is_bool)
    return valbool;
  else
    type_error(p_is_bool);
  return false; // dummy
}

//-----------------------------------------------------------------------
// MEMBER FUNCTIONS
//-----------------------------------------------------------------------

// FTHROW an error

void parameter_value::type_error(parameter_value_type expected) const {
  std::string text = "Parameter should be ";
  switch (expected) {
  case p_is_integer:
    text += "integer";
    break;

  case p_is_float:
    text += "floating point";
    break;

  case p_is_string:
    text += "string";
    break;

  case p_is_complex:
    text += "complex";
    break;

  case p_is_bool:
    text += "bool";
    break;

  default:
    Log.err("illegal parameter type in parameter_value::type_error");
  }

  text += " but is ";

  switch (expected) {
  case p_is_integer:
    text += "integer";
    break;

  case p_is_float:
    text += "floating point";
    break;

  case p_is_string:
    text += "string";
    break;

  case p_is_complex:
    text += "complex";
    break;

  case p_is_bool:
    text += "bool";
    break;

  case p_is_invalid:
    text += "invalid";

  default:
    Log.err("illegal parameter type in parameter_value::type_error");
  }

  text += ".\n";
  Log.err(text);
}

//=======================================================================
// write the parameter value into an ostream
//-----------------------------------------------------------------------

std::ostream &operator<<(std::ostream &out, const parameter_value &val) {
  switch (val.type()) {
  case p_is_integer:
    out << ((long)val);
    break;

  case p_is_float:
    out << ((double)val);
    break;

    //    case p_is_complex:
    //      out << ((complex&) val);
    //      break;

  case p_is_bool:
    out << ((bool)(val) ? "true" : "false");
    break;

  case p_is_string:
    out << '"' << ((const char *)val) << '"';
    break;

  case p_is_invalid:
    out << "invalid";
    break;

  default:
    Log.err("illegal parameter type in operator<<(ostream&,const "
	    "parameter_value&)");
  }

  return out;
}

} // namespace hydra
