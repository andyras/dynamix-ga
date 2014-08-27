#include "dynamixGAParams.hpp"

objectiveTypeFn dynamixGAParams::getObjectiveType() {
  if (gp.objectiveType.compare("single") == 0) {
    return singleObjective;
  }
  else if (gp.objectiveType.compare("double") == 0) {
    return doubleObjective;
  }
  else {
    std::cerr << "[" << __FUNCTION__ << "]: unrecognized objective type" << gp.objectiveType << std::endl;
    exit(-1);
  }
}

initializerFn dynamixGAParams::getInitializer() {
  // XXX: this is hacked a bit to just return what I hope to use in general, the
  // sensibleRandomInitializer.
  return sensibleRandomInitializer;
}

mutatorFn dynamixGAParams::getMutator() {
  return GA1DArraySensibleRandomMutator;
}