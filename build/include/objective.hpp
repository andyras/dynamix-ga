#ifndef __OBJECTIVE__
#define __OBJECTIVE__

#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>

#include <dynamix.hpp>
#include <ga-mpi/ga.h>
#include <propagate.hpp>

#include "gaparams.hpp"
#include "initializer.hpp"
#include "output.hpp"
#include "params.hpp"
#include "parser.hpp"

typedef float (*objectiveTypeFn)(GAGenome &); // function pointer to objective
typedef double (*objectiveFn)(Params *); // function pointer to subobjective

float singleObjective(GAGenome &);

float doubleObjective(GAGenome &);

double objAcceptorPeak(Params * p);

double objAcceptorAvg(Params * p);

double objAcceptorAvgAfterPeak(Params * p);

double objAcceptorFinal(Params * p);

objectiveTypeFn getObjectiveType(GAParams * p);

objectiveFn getObjective(GAParams * p);

#endif