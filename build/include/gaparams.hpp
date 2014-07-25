#ifndef __GAPARAMS__
#define __GAPARAMS__

#include <string>
#include <boost/mpi.hpp>

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
    int (*mutator)(GAGenome &, float); // function pointer to mutator

    // GA1DArrayGenome<double> bestGenome;

    bool firstEval = true;
    double bestScore = 0.0;
    std::vector<double> bestGenome;

    // upper and lower bounds for genes
    std::vector<double> lb;
    std::vector<double> ub;
  private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
      ar & objectiveType;
      ar & objective;
      ar & doubleObjective;
      ar & initializer;
      ar & minmax;

      ar & popsize;
      ar & pMut;
      ar & pCross;
      ar & convergence;

      ar & initializerFn;
      ar & objectiveFn;
      ar & mutator;

      ar & firstEval;
      ar & bestScore;
      ar & bestGenome;

      ar & lb;
      ar & ub;
      }
};

#endif