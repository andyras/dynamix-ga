#include "initializeGenome.hpp"

void sensibleRandomInitializer(GAGenome &g) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;
  GAParams * gp = (GAParams *) genome.userData();

  for (unsigned int ii = 0; ii < gp->paramsToChange.size(); ii++) {
    paramTuple ptc = gp->paramsToChange[ii];
    genome.gene(ii, GARandomFloat(std::get<1>(ptc), std::get<2>(ptc)));
  }

  return;
}