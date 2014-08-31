#pragma once

#include <tuple>
#include <string>
#include <boost/mpi.hpp>

#include <ga-mpi/ga.h>
#include <params.hpp>

typedef float (*objectiveTypeFn)(GAGenome &); // function pointer to objective
typedef double (*objectiveFn)(Params *); // function pointer to subobjective

typedef void (*initializerFn)(GAGenome &); // function pointer to initializer

typedef int (*mutatorFn)(GAGenome &, float); // function pointer to mutator

class GAParams {
  public:
    std::string objectiveType = "single";
    std::string doubleObjectiveType = "coherence";
    std::string initializer = "wavepacket";
    std::string minmax = "min";

    int popsize = 36;
    double pMut = 0.2;
    double pCross = 0.6;
    double convergence = 0.01;

    // a few system parameters not contained in Params
    double Da;  // spacing between acceptor states
    double Ea;  // acceptor lower band edge

    // variables below this line are not controlled by input file //////////////

    bool firstEval = true;
    double bestScore = 0.0;
    std::vector<double> bestGenome;

    // upper and lower bounds for genes
    std::vector<double> lb;
    std::vector<double> ub;

    std::vector< std::tuple<std::string, double, double> > paramsToChange;

    Params p;

  private:
    friend class boost::serialization::access;

    // this method needs to be updated whenever the data contents change.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
      ar & objectiveType;
      ar & doubleObjectiveType;
      ar & initializer;
      ar & minmax;

      ar & popsize;
      ar & pMut;
      ar & pCross;
      ar & convergence;

      ar & Da;
      ar & Ea;

      ar & firstEval;
      ar & bestScore;
      ar & bestGenome;

      ar & lb;
      ar & ub;

      ar & paramsToChange;

      ar & p;
    }
};