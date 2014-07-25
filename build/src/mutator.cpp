#include "mutator.hpp"

int GA1DArraySensibleRandomMutator (GAGenome &c, float pMut) {
  if (pMut <= 0.0) {
    return(0);
  }

  GA1DArrayGenome<double> &child=(GA1DArrayGenome<double> &)c;
  GAParams * gp = (GAParams *) child.userData();

  int nMut = 0;
  for (int ii = 0; ii < child.length(); ii++) {
    if (GAFlipCoin(pMut)) {
      child.gene(ii, GARandomFloat(gp->lb[ii], gp->ub[ii]));
      nMut++;
    }
  }
  return nMut;
}

mutatorFn getMutator(GAParams * p) {
  return GA1DArraySensibleRandomMutator;
}