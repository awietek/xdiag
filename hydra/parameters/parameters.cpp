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
// This file includes the class definitions for the parameter class
//=======================================================================

#include <algorithm>
#include <fstream>

#include "parameters.h"
#include "parameters_impl.h"
#include "parser.h"
#include "palm_exception.h"

using namespace std;
//=======================================================================
// Parameters
//
// an associative array of parameter values
//-----------------------------------------------------------------------
namespace hydra { namespace parameters {

ostream& operator<<(ostream& out, const Parameters& parms)
{
  for_each(parms.begin(),parms.end(),parameters_output(out));
  return out;
}

ostream& operator<<(ostream& out, const parameters_collection& parms_coll)
{
  for_each(parms_coll.begin(),parms_coll.end(),parameters_collection_output(out));
  return out;
}

istream& operator>>(istream& in, Parameters& parms)
{
  parser the_parser(in);
  the_parser >> parms;
  return in;
}

istream& operator>>(istream& in, parameters_collection& parms_coll)
{
  parser the_parser(in);
  the_parser >> parms_coll;
  return in;
}

// add new values by parsing the input

parser& operator>>(parser& in, Parameters& parms)
{
  char c;
  c=in.next_token_nows();
  do 
    {   
      // ignore extra semi-colons
      while(c==';')        
        c=in.next_token_nows();
      
      if(c==is_string)
        {       
          string s(in.value().get_string());
	  c=in.next_token_nows();
	  
	  if(c=='}') {
	    return in;
	  }

          if(c!='=')
            error ( parse_error("= expected in assignmanet while parsing Parameters") );
	  
          c=in.next_token_nows();
          switch(c)
            {
            case is_integer:
            case is_float:
            case is_string:
            case is_bool:
            case is_complex:
              parms[s] = in.value();
              break;
	      
            default:                            
              error(parse_error("invalid parameter value in input"));              
            }
	  
          // must be followed by a semicolon, comma or newline
          c=in.next_token();
          if((c!=';')&&(c!=',')&&(c!='\n'))
            error(parse_error("semicolon, comma or newline expected while parsing Parameters"));
          c=in.next_token_nows(); 
        }
      else 
        {
	  in.putback(c);
          return in;
        }
    } while (true);
}

parser& operator>>(parser& in, parameters_collection& parms_coll)
{
  char c;
  
  c=in.next_token_nows();
  do 
    {   
      // ignore extra semi-colons
      while(c==';')        
        c=in.next_token_nows();
      
      if(c==is_string)
        {       
          string s(in.value().get_string());
	  c=in.next_token_nows();
	  switch(c) {
	  case '{': 
	    // start parsing collection
	    {
	      Parameters parm;
	      in >> parm;
	      parms_coll[s]=parm;
	    }
	    c=in.next_token_nows();
	    break;
	  default:
	    error(parse_error("No global values allowed in collection"));
	  }
        }
      else 
        {
	  in.putback(c);
          return in;
        }
    } while (true);
}

Parameters::Parameters(parser& p)
{
  p >> (*this);
}

parameters_collection::parameters_collection(parser& p)
{
  p >> (*this);
}

Parameters read_parameters(std::string filename)
{
  // Open file and handle error
  std::ifstream File(filename.c_str());
  if(File.fail()) 
    {
      std::cerr << "Error in read_parameters: " 
		<< "Could not open file with filename ["
		<< filename << "] given. Abort." << std::endl;
      exit(EXIT_FAILURE);
    }

  parser prs(File);
  return Parameters(prs);
}
}} // namespace hydra::parameters
//=======================================================================
