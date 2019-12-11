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

#include "palm_exception.h"

//=======================================================================
// error
//
// function used to FTHROW an exception or to terminate if 
// exeptions are not implemented (unless the second argument is false)
//-----------------------------------------------------------------------
namespace hydra { namespace parameters {

void error(const std::exception& err,bool throw_it)
{
  throw err;
}

}} // namespace hydra::parameters
//=======================================================================
