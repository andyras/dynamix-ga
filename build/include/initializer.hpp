#ifndef __INITIALIZER__
#define __INITIALIZER__

#include <dynamix.hpp>
#include <ga-mpi/ga.h>

#include <assert.h>
#include "gaparams.hpp"
#include "params.hpp"

typedef void (*initializerFn)(GAGenome &); // function pointer to initializer

void gammasInitializer(GAGenome &g);

void init_gammas(GAGenome &c, Params * p);

void wavepacketInitializer(GAGenome &g);

void init_wavepacket(GAGenome &c, Params * p);

void wavepacketGammasInitializer(GAGenome &g);

void init_wavepacketGammas(GAGenome &c, Params * p);

void init_torsion(GAGenome &c, Params * p);

void sensibleRandomInitializer(GAGenome &g);

initializerFn getInitializer(GAParams * p);

#endif