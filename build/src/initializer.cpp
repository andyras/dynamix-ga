#include <dynamix.hpp>

#include "gaparams.hpp"
#include "initializer.hpp"

// #define DEBUG

void gammasInitializer(GAGenome &g) {

  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;

  genome.gene(0, GARandomFloat(0.001,0.05)); // gamma1
  genome.gene(1, GARandomFloat(0.001,0.05)); // gamma2
  genome.gene(2, GARandomFloat(0.001,0.05)); // gamma1_c

  return;
}

void init_gammas(GAGenome &c, Params * p) {
  // This function initializes the parameters which are to be changed in an
  // optimization over relaxation constants.

  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;

  // variables to be changed
  double g1 = genome.gene(0);
  double g2 = genome.gene(1);
  double g1_c = genome.gene(2);

  // assign parameters from GA /////////////////////////////////////////////////

  p->gamma1 = g1;
  p->gamma2 = g2;
  p->gamma1_c = g1_c;

  initialize(p);

  return;
}

void wavepacketInitializer(GAGenome &g) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;

  genome.gene(0, GARandomFloat(0.0,0.01)); // sigma
  genome.gene(1, GARandomFloat(0.0,0.01)); // mu

  return;
}

void init_wavepacket(GAGenome &c, Params * p) {
  // This function initializes the parameters which are to be changed in an
  // optimization over relaxation constants.

  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;

  // variables to be changed
  double bulkGaussSigma = genome.gene(0);
  double bulkGaussMu = genome.gene(1);

  // assign parameters from GA /////////////////////////////////////////////////

  p->bulkGaussSigma = bulkGaussSigma;
  p->bulkGaussMu = bulkGaussMu;

  initialize(p);

  return;
}

void wavepacketGammasInitializer(GAGenome &g) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;

  genome.gene(0, GARandomFloat(0.001, 0.05)); // gamma1
  genome.gene(1, GARandomFloat(0.001, 0.05)); // gamma2
  genome.gene(2, GARandomFloat(0.001, 0.05)); // gamma1_c
  genome.gene(3, GARandomFloat(0.0, 0.01)); // sigma
  genome.gene(4, GARandomFloat(0.0, 0.01)); // mu

  return;
}

void init_wavepacketGammas(GAGenome &c, Params * p) {
  // This function initializes the parameters which are to be changed in an
  // optimization over relaxation constants and wavepacket parameters.

  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;

  // variables to be changed
  double g1 = genome.gene(0);
  double g2 = genome.gene(1);
  double g1_c = genome.gene(2);
  double bulkGaussSigma = genome.gene(3);
  double bulkGaussMu = genome.gene(4);

  // assign parameters from GA /////////////////////////////////////////////////

  p->gamma1 = g1;
  p->gamma2 = g2;
  p->gamma1_c = g1_c;
  p->bulkGaussSigma = bulkGaussSigma;
  p->bulkGaussMu = bulkGaussMu;

  initialize(p);

  return;
}

void sensibleRandomInitializer(GAGenome &g) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;
  GAParams * gp = (GAParams *) genome.userData();

  for (unsigned int ii = 0; ii < genome.length(); ii++) {
    genome.gene(ii, GARandomFloat(gp->lb[ii], gp->ub[ii]));
  }

  return;
}

initializerFn getInitializer(GAParams * p) {
  return sensibleRandomInitializer;
}