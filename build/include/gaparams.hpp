#ifndef __GAPARAMS__
#define __GAPARAMS__

#include <string>

#include <params.hpp>

class GAParams : Params {
  public:
    std::string objectiveType = "single";
    std::string objective = "acceptorPeak";
    std::string doubleObjective = "coherence";
    std::string variables = "wavepacket";
    std::string minmax = "min";

    int popsize = 36;
    double pMut = 0.2;
    double pCross = 0.6;
    double convergence = 0.01;

    Params p;
};

#endif