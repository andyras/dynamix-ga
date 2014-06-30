#ifndef __INITIALIZER__
#define __INITIALIZER__

#include <ga-mpi/ga.h>

#include "params.hpp"

void init_gammas(GAGenome &c, Params * p);

void init_wavepacket(GAGenome &c, Params * p);

#endif