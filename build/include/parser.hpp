#ifndef __PARSER__
#define __PARSER__

#include <mpi.h>
#include <unistd.h>

#include "gaparams.hpp"

void assignGAParams(std::string inputFile, GAParams * p);

#endif
