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

extern GAParams gp; // global variable, ugh

float singleObjective(GAGenome &);

float doubleObjective(GAGenome &);

double objAcceptorPeak(Params * p);

double objAcceptorAvg(Params * p);

double objAcceptorAvgAfterPeak(Params * p);

double objAcceptorFinal(Params * p);

#endif
