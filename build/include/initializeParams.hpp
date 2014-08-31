#pragma once

#include <dynamix.hpp>
#include <ga-mpi/ga.h>
#include <params.hpp>

#include "gaparams.hpp"

void init_gammas(GAGenome &c, Params * p);

void init_wavepacket(GAGenome &c, Params * p);

void init_wavepacketGammas(GAGenome &c, Params * p);

void init_torsion(GAGenome &c, Params * p);

void initParamsToChange(GAGenome &g, Params * p);