#include "initializeGenome.hpp"

void sensibleRandomInitializer(GAGenome &g) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;
  GAParams * gp = (GAParams *) genome.userData();

  for (int ii = 0; ii < genome.length(); ii++) {
    genome.gene(ii, GARandomFloat(gp->lb[ii], gp->ub[ii]));
  }

  return;
}

initializerFn getInitializer(GAParams * p) {
  // XXX: this is hacked a bit to just return what I hope to use in general, the
  // sensibleRandomInitializer.
  return sensibleRandomInitializer;
}