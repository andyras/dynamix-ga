#ifndef __OBJTUPLE__
#define __OBJTUPLE__

#include <string>
#include <tuple>
#include <vector>

#include "objective.hpp"

class objTuple {
public:
  // internally the objTuple is two vectors, so that the object can be serialized
  std::vector<std::string> obs;
  std::vector<double> weights;

  void addObjWeight(std::string of, double w);
  std::vector< std::tuple<objectiveFn, double> > getObjWeights();

private:
  friend class boost::serialization::access;

  // this method needs to be updated whenever the data contents change.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & obs;
    ar & weights;
  }
};

#endif