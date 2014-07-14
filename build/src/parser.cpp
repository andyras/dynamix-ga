#include "parser.hpp"

void assignGAParams(std::string inputFile, GAParams * p) {
/* assigns GA params to the Params struct from the input file */
  std::string line;
  std::string input_param;
  std::string param_val;
  size_t equals_pos;
  size_t space_pos;

#ifdef DEBUG_INPUT_PARSER
  std::cout << "\nParsing file: " << inputFile << "\n";
#endif

  std::ifstream bash_in;  // declare input file stream

  bash_in.open(inputFile.c_str(), std::ios::in);  // open file as input stream
  if (bash_in.good() == false) {
    fprintf(stderr, "ERROR [Inputs]: file '%s' not available for reading\n", inputFile.c_str());
    exit(-1);
  }

  // read first line of input file
  getline (bash_in,line);

  // skip non-parameter lines
  while ( line != "## START INPUT PARAMETERS ##") {
    getline (bash_in,line);
  }

  while ( line != "## END INPUT PARAMETERS ##") {
    // skip comment lines
    if ( line.substr(0,1) == "#" ) {
      getline (bash_in,line);
      continue;
    }
    // find first equals sign
    equals_pos=line.find("=");
    // find first whitespace
    space_pos=(line.find(" ") > line.find("\t") ? line.find("\t") : line.find(" "));
    // parameter name is before equals sign
    input_param = line.substr(0,int(equals_pos));
    // parameter is after equals sign, before space
    param_val = line.substr(int(equals_pos)+1,int(space_pos)-int(equals_pos));
    // extract parameters
#ifdef DEBUG
    std::cout << "Parameter: " << input_param << std::endl << "New value: " << atof(param_val.c_str()) << std::endl;
#endif
    if (input_param == "objectiveType") { p->objectiveType = param_val; }
    else if (input_param == "objective") { p->objective = param_val; }
    else if (input_param == "doubleObjective") { p->doubleObjective = param_val; }
    else if (input_param == "variables") { p->variables = param_val; }
    else if (input_param == "minmax") { p->minmax = param_val; }
    else {  }
    getline (bash_in,line);
  }

  // close input file
  bash_in.close();

#ifdef DEBUG
  std::cout << std::endl;
  std::cout << "objectiveType is " << p->objectiveType << std::endl;
  std::cout << "objective is " << p->objective << std::endl;
  std::cout << "doubleObjective is " << p->doubleObjective << std::endl;
  std::cout << "variables is " << p->variables << std::endl;
  std::cout << "minmax is " << p->minmax << std::endl;
#endif

  // Error checking

  return;
}