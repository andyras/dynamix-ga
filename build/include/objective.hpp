#pragma once

#include <mpi.h>
#include <sys/types.h>
#include <tuple>
#include <unistd.h>

#include <dynamix.hpp>
#include <ga-mpi/ga.h>
#include <propagate.hpp>

#include "dynamixGAParams.hpp"
#include "initializeParams.hpp"
#include "output.hpp"
#include "params.hpp"
#include "subobjective.hpp"

float singleObjective(GAGenome &);

float doubleObjective(GAGenome &);