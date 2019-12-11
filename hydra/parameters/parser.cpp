//=======================================================================
//                 Alea library $Revision: 1.3 $
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
// This file includes functions used to parse the input parameters
//
// Changes needed:
// 
// It reads too far in case that other input should be performed
// on the stream 
//=======================================================================

#include <cctype>
#include <cstdio>
#include <cstdlib>

#include "parser.h"
#include "palm_exception.h"


//=======================================================================
// get the next (non whitespace) characzer
//-----------------------------------------------------------------------
namespace hydra { namespace parameters {

bool parser::isskip(char cc) 
{
  return isspace(cc)&&cc!='\n';
}

char parser::next_char() 
{ 
   do {
     c=in.get();
     } while (in&&isskip(c));
   return c;
}
  
char parser::eat_ws() 
{ 
  if(isspace(c)) in>>c; 
  return c;
}


//=======================================================================
// parse a std::string
//-----------------------------------------------------------------------

void parser::parse_string()
{
  // copy the read characters into the std::string until
  // either the closing quote, a new line, end of file or the maximum
  // length has been reached
  
  std::string the_string;
  do 
    {
      in.get(c);
      the_string += c;
    } while (c!='\n'&&c!='"'&&(!in.eof()));
  the_string.erase(the_string.length()-1,1); // remove last character
  
  // if there was a new line or end of file: error
  
  if(c!='"') 
    {
      error ( parse_error("unterminated string encountered while parsing") );
    }
    
  // prefetch the next character
  
  val=the_string;
  next_char();
}


//=======================================================================
// auxilliary function to parse an identifier
//-----------------------------------------------------------------------

parser::token parser::parse_ident()
{
   std::string the_string;
  the_string=c;

  // copy following alphanumeric characters or ' into the std::string
  do 
    {
      in.get(c);
      the_string+=c;
    } while ((isalnum(c)||c=='\'')&&(!in.eof()));
    
  the_string.erase(the_string.length()-1,1); // remove last character

  // special identifiers: yes, YES, true, TRUE
  
  if(the_string=="yes"||the_string=="true"||the_string=="YES"||the_string=="TRUE")
    {
      val = true;
      return is_bool;
    }

    
  // special identifiers: no, NO, false, FALSE
    
  if(the_string=="no"||the_string=="false"||the_string=="NO"||the_string=="FALSE")
    {
      val = false;
      return is_bool;
    }
  
  val = the_string;
  return is_string;
}
        

//=======================================================================
// auxilliary function to parse a number
//-----------------------------------------------------------------------

parser::token parser::parse_number()
{
  std::string the_string;

  // read an optional sign
  
  if(c=='+'||c=='-')
    {
      the_string+=c;
      c=in.get();
    }
  
  
  // read digits
  
  while(isdigit(c)&&(!in.eof()))
    {
      the_string+=c;
      c=in.get();
    }
    
  // decimal point or exponent: floating point number
    
  if(c=='.'||tolower(c)=='e')
    {
      // read fractional digits
      if(c=='.')
        {
          do 
            {
              the_string+=c;
              c=in.get();
            } while (isdigit(c)&&(!in.eof()));
        }
      
      // read exponential
      if(tolower(c)=='e')
        {
          the_string+=c;
          c=in.get();
          
          // optional sign
          if(c=='+'||c=='-')
            {
              the_string+=c;
              c=in.get();
            }

          // read digits of exponent
          while (isdigit(c)&&(!in.eof()))
            {
               the_string+=c;
               c=in.get();
            }
        }
        
      double the_number;
      // convert the std::string to a floating point number
      if(sscanf(the_string.c_str(),"%lf",&the_number)!=1)
        error ( parse_error("invalid number encountered while parsing") );

      val=the_number;
      
      // preread next character if necessary
      if(isskip(c)) next_char(); 
        
      return is_float;
    }
    
  // otherwise convert the std::string to an integer
  
  val=std::atol(the_string.c_str());

  // preread next character if necessary
  if(isskip(c))  next_char();
  
  return is_integer;
}

//=======================================================================
// auxilliary function to parse a complex number
//-----------------------------------------------------------------------

parser::token parser::parse_complex()
{
  next_char();
  parse_number();
  double re=val;
  if(c==',')
    {
      next_char();
      parse_number();
      double im=val;
      val = std::complex<double>(re,im);
    }
  else
    val=re;
  if(c!=')')
    error(parse_error("invalid complex number encountered while parsing") );
  return val.type();
}


//=======================================================================
// parse the next token
//-----------------------------------------------------------------------

// no whitespace token
char parser::next_token_nows() 
{
  char tok ;
  do {
    tok=next_token();
  } while (tok=='\n'||tok==' '||tok=='\t');
  return tok;
}

char parser::next_token()
{
  if(was_put_back)
    {
      was_put_back=false;
      return put_back_token;
    }

  if(in.eof() || c=='#') return eof; // special: EOF
  
  // starts with a letter: identifier
  if(isalpha(c)) 
    return parse_ident();

  // starts with a number, period, plus or minus: number
  if(isdigit(c)||c=='-'||c=='.'||c=='+') return parse_number();

  // starts with a ( : complex
  if(c=='(') return parse_complex();

  // starts with a ": std::string 
  if(c=='"') { parse_string(); return is_string;}
  
  // any other non white space character: special token
  char n=token(c);
  next_char();
  return n;
}

//=======================================================================
// put a token back
//-----------------------------------------------------------------------

bool parser::putback(char t)
{
  if(!was_put_back)
    {
      was_put_back=true;
      put_back_token=t;
      return true;
    }
   else
     return false;
}

}} // namespace hydra::parameters
