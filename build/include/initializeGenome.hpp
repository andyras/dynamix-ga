#ifndef __INITIALIZEGENOME__
#define __INITIALIZEGENOME__

// #include <dynamix.hpp>
#include <ga-mpi/ga.h>

#include "gaparams.hpp"

typedef void (*initializerFn)(GAGenome &); // function pointer to initializer

void sensibleRandomInitializer(GAGenome &g);

initializerFn getInitializer(GAParams * p);

#endif