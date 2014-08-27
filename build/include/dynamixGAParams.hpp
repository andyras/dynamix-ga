#ifndef __DYNAMIXGAPARAMS__
#define __DYNAMIXGAPARAMS__

#include <tuple>

#include "gaparams.hpp"
#include "initializeGenome.hpp"
#include "mutator.hpp"
#include "objTuple.hpp"
#include "objective.hpp"

class dynamixGAParams {
public:
  // data //////////////////////////////////////////////////////////////////////
  GAParams gp;
  objTuple ot;
  // end data //////////////////////////////////////////////////////////////////

  // methods ///////////////////////////////////////////////////////////////////
  objectiveTypeFn getObjectiveType();

  initializerFn getInitializer();

  mutatorFn getMutator();

private:
  friend class boost::serialization::access;

  // this method needs to be updated whenever the data contents change.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & gp;
    ar & ot;
  }
};

#endif