#ifndef __OBJECTIVE__
#define __OBJECTIVE__

#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>

#include <ga-mpi/ga.h>
#include <dynamix.hpp>
#include <propagate.hpp>

#include "params.hpp"
#include "gaparams.hpp"
#include "initializer.hpp"
#include "output.hpp"

float dynamixObjective(GAGenome &);

float dualObjective(GAGenome &);

double objAcceptorPeak(Params * p);

double objAcceptorAvg(Params * p);

double objAcceptorAvgAfterPeak(Params * p);

double objAcceptorFinal(Params * p);

#endif
