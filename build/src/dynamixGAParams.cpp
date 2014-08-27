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

objectiveFn dynamixGAParams::getObjective() {
  if (gp.objective.compare("acceptorPeak") == 0) {
    return objAcceptorPeak;
  }
  else if (gp.objective.compare("acceptorAvg") == 0) {
    return objAcceptorAvg;
  }
  else if (gp.objective.compare("acceptorAvgAfterPeak") == 0) {
    return objAcceptorAvgAfterPeak;
  }
  else if (gp.objective.compare("acceptorFinal") == 0) {
    return objAcceptorFinal;
  }
  else {
    std::cerr << "[" << __FUNCTION__ << "]: unrecognized objective" << gp.objective << std::endl;
    exit(-1);
  }
}

initializerFn dynamixGAParams::getInitializer() {
  // XXX: this is hacked a bit to just return what I hope to use in general, the
  // sensibleRandomInitializer.
  return sensibleRandomInitializer;
}