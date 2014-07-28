#include "parser.hpp"

// #define DEBUG
// #define DEBUGFAIL

void assignGAParams(std::string inputFile, GAParams * p) {
/* assigns GA params to the Params struct from the input file */
  std::string line;
  std::string input_param;
  std::string param_val;
  size_t equals_pos;
  size_t space_pos;

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

#ifdef DEBUG
  std::cout << "\nParsing file: " << inputFile << "\n";
#endif

  std::ifstream bash_in;  // declare input file stream

  bash_in.open(inputFile.c_str(), std::ios::in);  // open file as input stream
  unsigned int maxAttempts = 10;
  for (unsigned int ii = 1; ii <= maxAttempts; ii++) {
    if (bash_in.good() == false) {
      std::cerr << "ERROR [Inputs:" << mpi_rank << "]: file '" << inputFile <<
	"' not available for reading on attempt " << ii << "." << std::endl;
      if (bash_in.fail()) {
	std::cerr << "Problem is fail bit." << std::endl;
      }
      if (bash_in.bad()) {
	std::cerr << "Problem is fail bad." << std::endl;
      }
      if (bash_in.eof()) {
	std::cerr << "Problem is fail eof." << std::endl;
      }
      bash_in.close();
      usleep(1e6);
    }
    else if (ii == maxAttempts) {
      std::cerr << "dying after last attempt to read " << inputFile << std::endl;
      exit(-1);
    }
    else {
      break;
    }
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
    else if (input_param == "doubleObjectiveType") { p->doubleObjectiveType = param_val; }
    else if (input_param == "initializer") { p->initializer = param_val; }
    else if (input_param == "minmax") { p->minmax = param_val; }
    else if (input_param == "popsize") { p->popsize = atoi(param_val.c_str()); }
    else if (input_param == "pMut") { p->pMut = atof(param_val.c_str()); }
    else if (input_param == "pCross") { p->pCross = atof(param_val.c_str()); }
    else if (input_param == "convergence") { p->convergence = atof(param_val.c_str()); }
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

#ifdef DEBUGFAIL
  std::cerr << "[" << mpi_rank << "] fail bit for closing " << inputFile << " is " << bash_in.fail() << std::endl;
#endif

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
