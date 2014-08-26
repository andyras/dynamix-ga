#ifndef __DYNAMIXGAPARAMS__
#define __DYNAMIXGAPARAMS__

#include <tuple>

#include "gaparams.hpp"
#include "objective.hpp"
#include "objTuple.hpp"

class dynamixGAParams {
public:
  // data //////////////////////////////////////////////////////////////////////
  GAParams gp;
  objTuple ot;
  // end data //////////////////////////////////////////////////////////////////

  // methods ///////////////////////////////////////////////////////////////////
  objectiveTypeFn getObjectiveType();
  objectiveFn getObjective();

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