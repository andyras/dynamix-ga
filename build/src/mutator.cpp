#include "mutator.hpp"

int GA1DArraySensibleRandomMutator (GAGenome &c, float pMut) {
  if (pMut <= 0.0) {
    return(0);
  }

  GA1DArrayGenome<double> &child=(GA1DArrayGenome<double> &)c;
  GAParams * gp = (GAParams *) child.userData();

  int nMut = 0;
  double lb = 0.0;
  double ub = 0.0;

  for (unsigned int ii = 0; ii < gp->paramsToChange.size(); ii++) {
    if (GAFlipCoin(pMut)) {
      paramTuple ptc = gp->paramsToChange[ii];
      lb = std::get<1>(ptc);
      ub = std::get<2>(ptc);
      child.gene(ii, GARandomFloat(lb, ub));
      nMut++;
    }
  }

  return nMut;
}