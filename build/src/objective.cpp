#include "objective.hpp"

// #define DEBUG

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
  std::vector<double> tcprobs (p->numOutputSteps + 1, 0.0);

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

  f = 1.0 - findArrayMaximum(&(tcprobs[0]), tcprobs.size())/summ;
#ifdef DEBUG
  std::cout << "Value of f is 1.0 - " << findArrayMaximum(&(tcprobs[0]), p->numOutputSteps) << "/" << summ << " = " << f << std::endl <<
    "sigma: " << p->bulkGaussSigma << " mu: " << p->bulkGaussMu << " f: " << f << std::endl;
#endif

  return f;
}

double obj_Pcavg(Params * p) {
  // This function returns the average population on the QD states over the
  // course of the simulation

  double f = 0.0;

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

  double tcprob = 0.0;

  // at each time step, add populations in acceptor states
  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    if (p->wavefunction) {
      for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
        tcprob += pow(p->wfnt[ii*2*p->NEQ + jj],2) +
          pow(p->wfnt[ii*2*p->NEQ + jj + p->NEQ],2);
      }
    }
    else {
      for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
        tcprob += p->dmt[ii*2*p->NEQ2 + jj*p->NEQ + jj];
      }
    }
  }

  f = 1.0 - tcprob/summ/(p->numOutputSteps + 1);

#ifdef DEBUG
  std::cout << "Value of f is 1.0 - " << tcprob << "/" << summ << "/" << (p->numOutputSteps + 1) << " = " << f << std::endl <<
    "sigma: " << p->bulkGaussSigma << " mu: " << p->bulkGaussMu << " f: " << f << std::endl;
#endif

  return f;
}

double obj_Pcavg_after_peak(Params * p) {
  // This function returns 1.0 - the average population on the QD states over the
  // course of the simulation

  double f = 0.0;

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

  std::vector<double> tcprobs (p->numOutputSteps + 1, 0.0);

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

  unsigned int maxInd = findArrayMaximumIndex(&(tcprobs[0]), tcprobs.size());

  // add up populations on QD after peak
  double tcprob = 0.0;
  for (int ii = maxInd; ii <=p->numOutputSteps; ii++) {
    tcprob += tcprobs[ii];
  }

  f = 1.0 - tcprob/summ/(p->numOutputSteps + 1 - maxInd);

#ifdef DEBUG
  std::cout << "Value of f is 1.0 - " << tcprob << "/" << summ << "/" << (p->numOutputSteps + 1 - maxInd) << " = " << f << std::endl <<
    "sigma: " << p->bulkGaussSigma << " mu: " << p->bulkGaussMu << " f: " << f << std::endl;
#endif

  return f;
}

double obj_maxFinal(Params * p) {
  // This function returns 1.0 - the amount of population at the end of the
  // simulation. Notes:
  // 1. This sort of assumes a steady-state final population
  // 2. This is equivalent to minimizing BET and maximizing population, with
  // equal weights.

  double f = 0.0;

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

  double finalPop = 0.0;

  // at last time step, add populations in acceptor states
  int no = p->numOutputSteps;
  if (p->wavefunction) {
    for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
      finalPop += pow(p->wfnt[no*2*p->NEQ + jj],2) +
	pow(p->wfnt[no*2*p->NEQ + jj + p->NEQ],2);
    }
  }
  else {
    for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
      finalPop += p->dmt[no*2*p->NEQ2 + jj*p->NEQ + jj];
    }
  }

  f = 1.0 - finalPop;

#ifdef DEBUG
  std::cout << "sigma: " << p->bulkGaussSigma << " mu: " << p->bulkGaussMu << " f: " << f << std::endl;
#endif

  return f;
}
