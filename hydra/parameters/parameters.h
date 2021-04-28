//=======================================================================
//                 Alea library $Revision: 1.1.1.1 $
//          (c) 1994-1996 by Matthias Troyer and Beat Ammon
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
// This file includes the class definitions for the parameter classes
//=======================================================================

#include <complex>
#include <map>
#include <string>

#include <hydra/common.h>

namespace hydra {

class parser;
class parameter;

//=======================================================================
// types of parameter values
//-----------------------------------------------------------------------

typedef enum {
  p_is_invalid = 0,
  p_is_integer = 3,
  p_is_string = 4,
  p_is_float = 5,
  p_is_complex = 6,
  p_is_bool = 7,
  eof = 1
} parameter_value_type;

//=======================================================================
// parameter_value
//
// a parameter value, containing the variable name, type and value
//-----------------------------------------------------------------------

class parameter_value {
protected:
  std::string the_string; // the string value

  double valfloat;
  long valint;
  bool valbool;
  std::complex<double> valcomplex;

  parameter_value_type the_type; // the type

  void type_error(parameter_value_type) const;

public:
  // constructors: give name and optional value

  parameter_value();
  parameter_value(const parameter_value &);
  parameter_value(const std::string &);
  parameter_value(const long);
  parameter_value(const double);
  parameter_value(const bool);
  parameter_value(const std::complex<double> &);

  // assignment operators: assign a new value

  parameter_value &operator=(const parameter_value &);
  parameter_value &operator=(const char *);
  parameter_value &operator=(const std::string &);
  parameter_value &operator=(const long);
  parameter_value &operator=(const double);
  parameter_value &operator=(const bool);
  parameter_value &operator=(const std::complex<double> &);

  inline parameter_value &operator=(const int i) { return operator=((long)i); }

  inline parameter_value &operator=(const float x) {
    return operator=((double)x);
  }

  // conversion operators

  operator const std::string &() const;
  operator const char *() const;
  operator long() const;
  operator double() const;
  operator int() const;
  operator float() const;
  operator bool() const;
  operator std::complex<double>() const;

  const std::string &get_string() const {
    return operator const std::string &();
  }

  long get_integer() const { return operator long(); }

  double get_number() const { return operator double(); }

  std::complex<double> get_complex() const {
    return operator std::complex<double>();
  }

  bool get_boolean() const { return operator bool(); }

  // return the type
  parameter_value_type type() const { return the_type; }
  int is_numeric() const {
    return (type() == p_is_integer) || (type() == p_is_float) ||
           (type() == p_is_bool) || (type() == p_is_complex);
  }
};

//=======================================================================
// Parameters
//
// an associative array of parameter values
//-----------------------------------------------------------------------

class Parameters : public std::map<std::string, parameter_value> {
public:
  Parameters() : std::map<std::string, parameter_value>(){};
  Parameters(const Parameters &p) : std::map<std::string, parameter_value>(p){};
  Parameters(parser &);
  // copy/replace Parameters into this one
  Parameters &operator<<(const Parameters &p) {
    insert(p.begin(), p.end());
    return *this;
  }

  // return type of parameter with given name
  inline bool defined(const std::string &n) const {
    return (this->find(n) != this->end());
  }
  inline bool undefined(const std::string &n) const {
    return (this->find(n) == this->end());
  }
};

class parameters_collection : public std::map<std::string, Parameters> {
public:
  parameters_collection() : std::map<std::string, Parameters>(){};
  parameters_collection(const parameters_collection &pc)
      : std::map<std::string, Parameters>(pc){};
  parameters_collection(parser &);

  // copy/replace parameters into this one
  parameters_collection &operator<<(const parameters_collection &pc) {
    insert(pc.begin(), pc.end());
    return *this;
  }

  // return type of parameter with given name
  inline bool defined(const std::string &n) {
    return (this->find(n) != this->end());
  }
  inline bool undefined(const std::string &n) {
    return (this->find(n) == this->end());
  }
};

//=======================================================================
// functions for reading parameters from a parser
//-----------------------------------------------------------------------

extern parser &operator>>(parser &, Parameters &);
extern parser &operator>>(parser &, parameters_collection &);

extern std::ostream &operator<<(std::ostream &, const parameters_collection &);
extern std::ostream &operator<<(std::ostream &, const Parameters &);
extern std::ostream &operator<<(std::ostream &, const parameter_value &);

Parameters read_parameters(std::string filename);

} // namespace hydra
