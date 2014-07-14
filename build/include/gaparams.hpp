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

    Params p;
};

#endif