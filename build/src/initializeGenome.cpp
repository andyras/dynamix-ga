#include "initializeGenome.hpp"

void sensibleRandomInitializer(GAGenome &g) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;
  GAParams * gp = (GAParams *) genome.userData();

  for (int ii = 0; ii < genome.length(); ii++) {
    genome.gene(ii, GARandomFloat(gp->lb[ii], gp->ub[ii]));
  }

  return;
}