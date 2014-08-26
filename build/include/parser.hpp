#ifndef __PARSER__
#define __PARSER__

#include <boost/algorithm/string.hpp>
#include <mpi.h>
#include <unistd.h>

#include "dynamixGAParams.hpp"

void assignGAParams(std::string inputFile, dynamixGAParams * dgp);

#endif