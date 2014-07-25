#ifndef __MUTATOR__
#define __MUTATOR__

#include <ga-mpi/ga.h>

#include "gaparams.hpp"

typedef int (*mutatorFn)(GAGenome &, float); // function pointer to mutator

int GA1DArraySensibleRandomMutator (GAGenome &c, float pMut);

mutatorFn getMutator(GAParams * p);

#endif