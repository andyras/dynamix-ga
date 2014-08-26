#include "objTuple.hpp"

// adds an objective function and weight to the corresponding lists
void objTuple::addObjWeight(std::string of, double w) {
  obs.emplace_back(of);
  weights.emplace_back(w);

  return;
}

// returns a vector of tuples with objective function pointers and weights
std::vector< std::tuple<objectiveFn, double> > objTuple::getObjWeights() {
  std::vector< std::tuple<objectiveFn, double> > objectives;

  // first token determines objective function pointer
  for (unsigned int ii = 0; ii < obs.size(); ii++) {
    if (obs[ii].compare("acceptorPeak") == 0) {
      objectives.emplace_back(std::make_tuple(objAcceptorPeak, weights[ii]));
    }
    else if (obs[ii].compare("acceptorAvg") == 0) {
      objectives.emplace_back(std::make_tuple(objAcceptorAvg, weights[ii]));
    }
    else if (obs[ii].compare("acceptorAvgAfterPeak") == 0) {
      objectives.emplace_back(std::make_tuple(objAcceptorAvgAfterPeak, weights[ii]));
    }
    else if (obs[ii].compare("acceptorFinal") == 0) {
      objectives.emplace_back(std::make_tuple(objAcceptorFinal, weights[ii]));
    }
    else {
      std::cerr << "ERROR [" << __FUNCTION__ << "]: objective function " << obs[ii] << " not recognized." << std::endl;
    }
  }

  return objectives;
}