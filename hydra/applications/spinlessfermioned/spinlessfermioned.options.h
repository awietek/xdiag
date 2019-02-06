#ifndef HYDRA_APPLICATIONS_SPINLESSFERMIONED_OPTIONS_H_
#define HYDRA_APPLICATIONS_SPINLESSFERMIONED_OPTIONS_H_

#include <string>
#include "clara.hpp"

void parse_cmdline(double& t, double& V, std::string& outfile, std::string& latticefile, 
		   int& np, std::string& representation, int& argc, char** argv)
{
  
  bool showhelp = false;
  auto parser =
    clara::Opt(t, "t")["-t"]["--t"]("hopping constant t") |
    clara::Opt(V, "V")["-V"]["--V"]("interaction strength V") |
    clara::Opt(outfile, "outfile")["-o"]["--outfile"]("name of outfile") |
    clara::Opt(latticefile, "latticefile")["-l"]["--latticefile"]("name of latticefile") |
    clara::Opt(np, "np")["-p"]["--np"]("number of spinless fermions (optional, default: half filled)") |
    clara::Opt(representation, "representation")["-r"]["--representation"]("space group representation (optional: default no space group symmetries applied)") |
    clara::Help(showhelp);


  auto cmd_args = parser.parse(clara::Args(argc,argv));
  if( !cmd_args ) 
    {
      std::cerr << "Error in command line: " << cmd_args.errorMessage() << std::endl;
      exit(EXIT_FAILURE);
    }
  else if (showhelp) 
    {
      parser.writeToStream(std::cout);
      exit(EXIT_SUCCESS);
    }
  else
    {
      std::cout <<
	"t             : " << t << std::endl <<
	"V             : " << V << std::endl <<
	"outfile       : " << outfile << std::endl <<
	"latticefile   : " << latticefile << std::endl <<
	"np            : " << np << std::endl <<
	"representation: " << representation << std::endl <<
	"-----" << std::endl;
    }
}


#endif
