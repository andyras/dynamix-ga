#ifndef __GAPARAMS__
#define __GAPARAMS__

#include <string>

#include <ga-mpi/ga.h>

#include "initializer.hpp"

class GAParams {
  public:
    std::string objectiveType = "single";
    std::string objective = "acceptorPeak";
    std::string doubleObjective = "coherence";
    std::string initializer = "wavepacket";
    std::string minmax = "min";

    int popsize = 36;
    double pMut = 0.2;
    double pCross = 0.6;
    double convergence = 0.01;

    // variables below this line are not controlled by input file //////////////

    void (*initializerFn)(GAGenome &); // function pointer to initializer
    float (*objectiveFn)(GAGenome &); // function pointer to objective

    // GA1DArrayGenome<double> bestGenome;

    bool firstEval = true;
    double bestScore = 0.0;
    std::vector<double> bestGenome;
};

#endif