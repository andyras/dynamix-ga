#pragma once

#include <boost/algorithm/string.hpp>
#include <mpi.h>
#include <unistd.h>

#include "dynamixGAParams.hpp"
#include "serialize_tuple.hpp"

void assignGAParams(std::string inputFile, dynamixGAParams * dgp);