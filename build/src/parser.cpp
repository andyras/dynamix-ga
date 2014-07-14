#include "parser.hpp"

// #define DEBUG

void assignGAParams(std::string inputFile, GAParams * p) {
/* assigns GA params to the Params struct from the input file */
  std::string line;
  std::string input_param;
  std::string param_val;
  size_t equals_pos;
  size_t space_pos;

#ifdef DEBUG
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
    if (input_param == "objectiveType") { p->objectiveType = param_val; }
    else if (input_param == "objective") { p->objective = param_val; }
    else if (input_param == "doubleObjective") { p->doubleObjective = param_val; }
    else if (input_param == "variables") { p->variables = param_val; }
    else if (input_param == "minmax") { p->minmax = param_val; }
    else if (input_param == "popsize") { p->popsize = stoi(param_val); }
    else if (input_param == "pMut") { p->pMut = stod(param_val); }
    else if (input_param == "pCross") { p->pCross = stod(param_val); }
    else if (input_param == "convergence") { p->convergence = stod(param_val); }
    else {
      std::cerr << "WARNING: unknown input parameter " << input_param << std::endl;
    }
#ifdef DEBUG
    std::cout << "Parameter: " << input_param << std::endl << "New value: " << param_val << std::endl;
#endif
    getline (bash_in,line);
  }

  // close input file
  bash_in.close();

  // Error checking

  // popsize / mpi_tasks must be an integer
  int mpi_tasks = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
  int newPopSize = mpi_tasks * int((double)p->popsize/(double)mpi_tasks+0.99999);
  if (p->popsize % mpi_tasks != 0) {
    std::cerr << "WARNING: MPI task number (" << mpi_tasks <<
      ") does not divide evenly into population size (" << p->popsize <<
      "). Setting population size to " << newPopSize << std::endl;
    p->popsize = newPopSize;
  }

  return;
}