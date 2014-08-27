#ifndef __OBJECTIVE__
#define __OBJECTIVE__

#include <mpi.h>
#include <sys/types.h>
#include <tuple>
#include <unistd.h>

#include <dynamix.hpp>
#include <ga-mpi/ga.h>
#include <propagate.hpp>

#include "dynamixGAParams.hpp"
#include "initializer.hpp"
#include "output.hpp"
#include "params.hpp"
#include "subobjective.hpp"

float singleObjective(GAGenome &);

float doubleObjective(GAGenome &);

#endif