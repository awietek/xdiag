//=======================================================================
//                 Osiris library $Revision: 1.1.1.1 $
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

#ifndef ___OSIRIS_EXCEPTION___
#define ___OSIRIS_EXCEPTION___

//=======================================================================
// This file includes the various types of std::exceptions that can be thrown
// by the Osiris library.
//=======================================================================

//#include <ios>
#include <stdexcept>
#include <exception>
#include <string>

namespace hydra { namespace parameters {
//=======================================================================
// error
//
// function used to throw an std::exception or to terminate if 
// exeptions are not implemented (unless the second argument is false)
//-----------------------------------------------------------------------
   
void error(const std::exception& err,bool=true);
  

//=======================================================================
// Assert
//
// templates to assert a condition and throw an error of not valid
// taken from Stroustrup
//-----------------------------------------------------------------------

template <class X, class A>
inline void Assert(A assertion)
{
  if(!assertion) error(X());
}


template <class A, class E>
inline void Assert(A assertion, E except)
{
  if(!assertion) error(except);
}


//=======================================================================
// comm_error
//
// general std::exception in communication
//-----------------------------------------------------------------------

class comm_error : public std::runtime_error
{
public:
  comm_error(const std::string& s) throw() : runtime_error(s) {};
};


//=======================================================================
// os_error
//
// std::exception while calling low level OS function
//-----------------------------------------------------------------------

class os_error : public std::runtime_error
{
public:
  os_error(const std::string& s) throw() : runtime_error(s) {};
};


//=======================================================================
// invariants_error
//
// an invariant was broken
//-----------------------------------------------------------------------

class invariants_error : public std::logic_error
{
public:
  invariants_error(const std::string& s) throw() : logic_error(s) {};
};


//=======================================================================
// zero_pointer_error
//
// a pointer was zero, but it should point to a vaild object
//-----------------------------------------------------------------------

class zero_pointer_error : public std::runtime_error
{
public:
  zero_pointer_error(const std::string& s) throw() : runtime_error(s) {};
};


//=======================================================================
// default_error
//
// A branch in a switch or if statement was reached that should never
// be reached
//-----------------------------------------------------------------------

class default_error : public std::out_of_range
{
public:
  default_error(const std::string& s) throw() : out_of_range(s) {};
};


//=======================================================================
// illegal_code
//
// this code should never be executed
//-----------------------------------------------------------------------

class illegal_code : public std::logic_error
{
public:
  illegal_code(const std::string& s) throw() : logic_error(s) {};
};


//=======================================================================
// not_implemented_error
//
// this functionality has not been implemented yet
//-----------------------------------------------------------------------

class not_implemented : public std::logic_error
{
public:
  not_implemented(const std::string& s) throw() : logic_error(s) {};
};


//=======================================================================
// dump_error
//
// an error occured while reading from or writing to a dump
//-----------------------------------------------------------------------

class dump_error : public std::runtime_error
{
public:
  dump_error(const std::string& s) throw() : runtime_error(s) {};
};


//=======================================================================
// file_open_error
//
//-----------------------------------------------------------------------

/*
class file_open_error : public ios_base::failure
{
public:
  file_open_error(const std::string& s) throw() : ios_base::failure(s) {};
};
*/


//=======================================================================
// version_error
//
// cannot read this version of a class
//-----------------------------------------------------------------------

class version_error : public std::range_error
{
public:
  version_error(const std::string& s) throw() : range_error(s) {};
};


//=======================================================================
// dump_type_error
//
// wrong type of dump
//-----------------------------------------------------------------------

class dump_type_error : public std::range_error
{
public:
  dump_type_error(const std::string& s) throw() : range_error(s) {};
};


//=======================================================================
// parameter_type_error
//
// wrong type of parameter
//-----------------------------------------------------------------------

class parameter_type_error : public std::runtime_error
{
public:
  parameter_type_error(const std::string& s) throw() : runtime_error(s) {};
};


//=======================================================================
// parse_error
//
// error occured whiele parsing
//-----------------------------------------------------------------------

class parse_error : public std::runtime_error
{
public:
  parse_error(const std::string& s) throw() : runtime_error(s) {};
};

//=======================================================================
// math_error
//
// error occured in mathematical calculation
//-----------------------------------------------------------------------

class math_error : public std::runtime_error
{
public:
  math_error(const std::string& s) throw() : runtime_error(s) {};
};


}} // namespace hydra::parameters

//=======================================================================
#endif
//=======================================================================
