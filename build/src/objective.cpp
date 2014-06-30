#include "objective.hpp"

#define DEBUG

double obj_tcpeak(Params * p) {
  // This function returns the peak total population on the QD states
  double f = 0.0; // objective function return value

  // get initial population (in case it is > 1)
  double summ = 0.0;
  if (p->wavefunction) {
    for (int ii = 0; ii < p->NEQ; ii++) {
      summ += pow(p->wfnt[ii],2) + pow(p->wfnt[ii+p->NEQ],2);
    }
  }
  else {
    for (int ii = 0; ii < p->NEQ2; ii += (p->NEQ+1)) {
      summ += p->dmt[ii];
    }
  }

  // create array of populations
  std::vector<double> tcprobs (p->numOutputSteps, 0.0);

  // at each time step, add populations in acceptor states
  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    if (p->wavefunction) {
      for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
        tcprobs[ii] += pow(p->wfnt[ii*2*p->NEQ + jj],2) +
          pow(p->wfnt[ii*2*p->NEQ + jj + p->NEQ],2);
      }
    }
    else {
      for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
        tcprobs[ii] += p->dmt[ii*2*p->NEQ2 + jj*p->NEQ + jj];
      }
    }
  }

  f = findArrayMaximum(&(tcprobs[0]), p->Nc)/summ;
#ifdef DEBUG
  std::cout << "Value of f is " << f << std::endl;
#endif

  return (1.0 - f);
}