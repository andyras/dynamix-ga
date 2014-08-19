#include "objective.hpp"

// #define DEBUG

float singleObjective(GAGenome &c) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;
  GAParams * gp = (GAParams *) genome.userData();

#ifdef DEBUG
  print1DGenes(c);
#endif

  // get MPI rank and size
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int pid = getpid();

  // output value
  float output = 1.0;

  // Struct of parameters //////////////////////////////////////////////////////

  Params p = gp->p;

  // set number of processors for OpenMP ///////////////////////////////////////

  omp_set_num_threads(p.nproc);
  mkl_set_num_threads(p.nproc);

  // assign GA parameters
  // this function is independent of the objective you are using. It determines
  // what the relevant parameters are for the optimization.
  if (gp->initializer.compare("g1g2g1_c") == 0) {
    init_gammas(c, &p);
  }
  else if (gp->initializer.compare("wavepacket") == 0) {
    init_wavepacket(c, &p);
  }
  else if (gp->initializer.compare("wavepacketGammas") == 0) {
    init_wavepacketGammas(c, &p);
  }
  else if (gp->initializer.compare("torsion") == 0) {
    init_torsion(c, &p);
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "variable set" <<
      gp->initializer << "not recognized." << std::endl;
    exit(-1);
  }

  makePlots(&p);
  propagate(&p);

  // calculate value of objective function /////////////////////////////////////
#ifdef DEBUG
  std::cout << "[" << pid << ":" << rank << "] " <<
    "Calculating value of objective function..." << std::endl;
#endif
  // set 'output' according to an objective function ///////////////////////////
  if (gp->objective.compare("acceptorPeak") == 0) {
    output = objAcceptorPeak(&p);
  }
  else if (gp->objective.compare("acceptorAvg") == 0) {
    output = objAcceptorAvg(&p);
  }
  else if (gp->objective.compare("acceptorAvgAfterPeak") == 0) {
    output = objAcceptorAvgAfterPeak(&p);
  }
  else if (gp->objective.compare("acceptorFinal") == 0) {
    output = objAcceptorFinal(&p);
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "objective function" <<
      gp->objective << "not recognized." << std::endl;
    exit(-1);
  }

  std::cout << "[" << pid << ":" << rank << "] " << "Objective: " << output << " Genome:";
  for (unsigned int ii = 0; ii < (unsigned int) genome.length(); ii++) {
    std::cout << " " << genome.gene(ii);
  }
  std::cout << std::endl;

  // initialize best genome if this is the first call to the objective function
  if (gp->firstEval && (rank == 0)) {
    gp->firstEval = false;
#ifdef DEBUG
    std::cout << "[" << pid << ":" << rank << "] " << "first call to objective" << std::endl;
#endif
  }
  else {
#ifdef DEBUG
    std::cout << "[" << pid << ":" << rank << "] " << "not first call to objective" << std::endl;
#endif
  }

  // update best genome if this call is the best ///////////////////////////////
  bool isBest = false;
  if (((gp->minmax.compare("min") == 0) && (output < gp->bestScore)) ||
       ((gp->minmax.compare("max") == 0) && (output > gp->bestScore))) {
      isBest = true;
  }

  if (isBest) {
    gp->bestScore = output;
    for (unsigned int ii = 0; ii < gp->bestGenome.size(); ii++) {
      gp->bestGenome[ii] = genome.gene(ii);
    }
#ifdef DEBUG
    std::cout << "[" << pid << ":" << rank << "] " << "New best score: " << output << " Genome:";
    for (unsigned int ii = 0; ii < gp->bestGenome.size(); ii++) {
      std::cout << " " << gp->bestGenome[ii];
    }
    std::cout << std::endl;
#endif
  }
#ifdef DEBUG
  else {
    std::cout << "[" << pid << ":" << rank << "] " << "Not best score." << std::endl;
  }
#endif

  return output;
}

float doubleObjective(GAGenome &c) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;
  GAParams * gp = (GAParams *) genome.userData();

#ifdef DEBUG
  print1DGenes(c);
#endif

  int pid = getpid();
  // get MPI rank and size
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

#ifdef DEBUG
  std::cout << "rank: " << rank << " of " << size << "\n";
#endif

  // output value
  float output1 = 1.0;  // objective value 1
  float output2 = 1.0;  // objective value 2

  // Struct of parameters //////////////////////////////////////////////////////

  Params p = gp->p;

  // set number of processors for OpenMP ///////////////////////////////////////

  omp_set_num_threads(p.nproc);
  mkl_set_num_threads(p.nproc);

  // assign GA parameters
  // this function is independent of the objective you are using. It determines
  // what the relevant parameters are for the optimization.
  if (gp->initializer.compare("g1g2g1_c") == 0) {
    init_gammas(c, &p);
  }
  else if (gp->initializer.compare("wavepacket") == 0) {
    init_wavepacket(c, &p);
  }
  else if (gp->initializer.compare("wavepacketGammas") == 0) {
    init_wavepacketGammas(c, &p);
  }
  else if (gp->initializer.compare("torsion") == 0) {
    init_torsion(c, &p);
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "variable set" <<
      gp->initializer << "not recognized." << std::endl;
    exit(-1);
  }

  // coherent propagation //////////////////////////////////////////////////////

  p.coherent = 1;
  p.outputDir = "./coh/";
  makePlots(&p);
  initWavefunction(&p);
  propagate(&p);

  // calculate value of objective function /////////////////////////////////////
#ifdef DEBUG
  std::cout << "[" << pid << ":" << rank << "] " <<
    "Calculating value of objective function..." << std::endl;
#endif
  if (gp->objective.compare("acceptorPeak") == 0) {
    output1 = objAcceptorPeak(&p);
  }
  else if (gp->objective.compare("acceptorAvg") == 0) {
    output1 = objAcceptorAvg(&p);
  }
  else if (gp->objective.compare("acceptorAvgAfterPeak") == 0) {
    output1 = objAcceptorAvgAfterPeak(&p);
  }
  else if (gp->objective.compare("acceptorFinal") == 0) {
    output1 = objAcceptorFinal(&p);
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "objective function" <<
      gp->objective << "not recognized." << std::endl;
    exit(-1);
  }

  std::cout << "[" << pid << ":" << rank << "] " << "Objective 1: " << output1 << " Genome:";
  for (unsigned int ii = 0; ii < (unsigned int) genome.length(); ii++) {
    std::cout << " " << genome.gene(ii);
  }
  std::cout << std::endl;

  // incoherent propagation ////////////////////////////////////////////////////
  p.coherent = 0;
  p.outputDir = "./inc/";
  makePlots(&p);
  initWavefunction(&p);
  propagate(&p);

  // calculate value of objective function /////////////////////////////////////
#ifdef DEBUG
  std::cout << "[" << pid << ":" << rank << "] " <<
    "Calculating value of objective function again..." << std::endl;
#endif
  if (gp->objective.compare("acceptorPeak") == 0) {
    output2 = objAcceptorPeak(&p);
  }
  else if (gp->objective.compare("acceptorAvg") == 0) {
    output2 = objAcceptorAvg(&p);
  }
  else if (gp->objective.compare("acceptorAvgAfterPeak") == 0) {
    output2 = objAcceptorAvgAfterPeak(&p);
  }
  else if (gp->objective.compare("acceptorFinal") == 0) {
    output2 = objAcceptorFinal(&p);
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "objective function" <<
      gp->objective << "not recognized." << std::endl;
    exit(-1);
  }

  std::cout << "[" << pid << ":" << rank << "] " << "Objective 2: " << output2 << " Genome:";
  for (unsigned int ii = 0; ii < (unsigned int) genome.length(); ii++) {
    std::cout << " " << genome.gene(ii);
  }
  std::cout << std::endl;

  double output = fabs(output1 - output2);
  std::cout << "[" << pid << ":" << rank << "] " << "Combined objective: " << output << std::endl;

  // update best genome if this call is the best ///////////////////////////////
  bool isBest = false;
  if (((gp->minmax.compare("min") == 0) && (output < gp->bestScore)) ||
       ((gp->minmax.compare("max") == 0) && (output > gp->bestScore))) {
      isBest = true;
  }

  if (isBest) {
    gp->bestScore = output;
    for (unsigned int ii = 0; ii < gp->bestGenome.size(); ii++) {
      gp->bestGenome[ii] = genome.gene(ii);
    }
#ifdef DEBUG
    std::cout << "[" << pid << ":" << rank << "] " << "New best score: " << output << " Genome:";
    for (unsigned int ii = 0; ii < gp->bestGenome.size(); ii++) {
      std::cout << " " << gp->bestGenome[ii];
    }
    std::cout << std::endl;
#endif
  }
#ifdef DEBUG
  else {
    std::cout << "[" << pid << ":" << rank << "] " << "Not best score." << std::endl;
  }
#endif

  return output;
}

double objAcceptorPeak(Params * p) {
  // This function returns the peak total population on the QD states
  double f = 0.0; // objective function return value

  // get initial population (in case it is > 1)
  double summ = 0.0;
  if (p->wavefunction) {
    for (int ii = 0; ii < p->NEQ; ii++) {
      summ += pow(p->wfnt[ii],2) + pow(p->wfnt[ii+p->NEQ],2);
    }
  }
  else {
    for (int ii = 0; ii < p->NEQ2; ii += (p->NEQ+1)) {
      summ += p->dmt[ii];
    }
  }

  // create array of populations
  std::vector<double> tcprobs (p->numOutputSteps + 1, 0.0);

  // at each time step, add populations in acceptor states
  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    if (p->wavefunction) {
      for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
        tcprobs[ii] += pow(p->wfnt[ii*2*p->NEQ + jj],2) +
          pow(p->wfnt[ii*2*p->NEQ + jj + p->NEQ],2);
      }
    }
    else {
      for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
        tcprobs[ii] += p->dmt[ii*2*p->NEQ2 + jj*p->NEQ + jj];
      }
    }
  }

  f = findArrayMaximum(&(tcprobs[0]), tcprobs.size())/summ;

  return f;
}

double objAcceptorAvg(Params * p) {
  // This function returns the average population on the QD states over the
  // course of the simulation

  double f = 0.0;

  // get initial population (in case it is > 1)
  double summ = 0.0;
  if (p->wavefunction) {
    for (int ii = 0; ii < p->NEQ; ii++) {
      summ += pow(p->wfnt[ii],2) + pow(p->wfnt[ii+p->NEQ],2);
    }
  }
  else {
    for (int ii = 0; ii < p->NEQ2; ii += (p->NEQ+1)) {
      summ += p->dmt[ii];
    }
  }

  double tcprob = 0.0;

  // at each time step, add populations in acceptor states
  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    if (p->wavefunction) {
      for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
        tcprob += pow(p->wfnt[ii*2*p->NEQ + jj],2) +
          pow(p->wfnt[ii*2*p->NEQ + jj + p->NEQ],2);
      }
    }
    else {
      for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
        tcprob += p->dmt[ii*2*p->NEQ2 + jj*p->NEQ + jj];
      }
    }
  }

  f = tcprob/summ/(p->numOutputSteps + 1);

  return f;
}

double objAcceptorAvgAfterPeak(Params * p) {
  // This function returns 1.0 - the average population on the QD states over the
  // course of the simulation

  double f = 0.0;

  // get initial population (in case it is > 1)
  double summ = 0.0;
  if (p->wavefunction) {
    for (int ii = 0; ii < p->NEQ; ii++) {
      summ += pow(p->wfnt[ii],2) + pow(p->wfnt[ii+p->NEQ],2);
    }
  }
  else {
    for (int ii = 0; ii < p->NEQ2; ii += (p->NEQ+1)) {
      summ += p->dmt[ii];
    }
  }

  std::vector<double> tcprobs (p->numOutputSteps + 1, 0.0);

  // at each time step, add populations in acceptor states
  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    if (p->wavefunction) {
      for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
        tcprobs[ii] += pow(p->wfnt[ii*2*p->NEQ + jj],2) +
          pow(p->wfnt[ii*2*p->NEQ + jj + p->NEQ],2);
      }
    }
    else {
      for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
        tcprobs[ii] += p->dmt[ii*2*p->NEQ2 + jj*p->NEQ + jj];
      }
    }
  }

  unsigned int maxInd = findArrayMaximumIndex(&(tcprobs[0]), tcprobs.size());

  // add up populations on QD after peak
  double tcprob = 0.0;
  for (int ii = maxInd; ii <=p->numOutputSteps; ii++) {
    tcprob += tcprobs[ii];
  }

  f = tcprob/summ/(p->numOutputSteps + 1 - maxInd);

  return f;
}

double objAcceptorFinal(Params * p) {
  // This function returns 1.0 - the amount of population at the end of the
  // simulation. Notes:
  // 1. This sort of assumes a steady-state final population
  // 2. This is equivalent to minimizing BET and maximizing population, with
  // equal weights.

  double f = 0.0;

  // get initial population (in case it is > 1)
  double summ = 0.0;
  if (p->wavefunction) {
    for (int ii = 0; ii < p->NEQ; ii++) {
      summ += pow(p->wfnt[ii],2) + pow(p->wfnt[ii+p->NEQ],2);
    }
  }
  else {
    for (int ii = 0; ii < p->NEQ2; ii += (p->NEQ+1)) {
      summ += p->dmt[ii];
    }
  }

  double finalPop = 0.0;

  // at last time step, add populations in acceptor states
  int no = p->numOutputSteps;
  if (p->wavefunction) {
    for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
      finalPop += pow(p->wfnt[no*2*p->NEQ + jj],2) +
	pow(p->wfnt[no*2*p->NEQ + jj + p->NEQ],2);
    }
  }
  else {
    for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
      finalPop += p->dmt[no*2*p->NEQ2 + jj*p->NEQ + jj];
    }
  }

  f = finalPop;

  return f;
}


objectiveTypeFn getObjectiveType(GAParams * p) {
  if (p->objectiveType.compare("single") == 0) {
    return singleObjective;
  }
  else if (p->objectiveType.compare("double") == 0) {
    return doubleObjective;
  }
  else {
    std::cerr << "[" << __FUNCTION__ << "]: unrecognized objective type" << p->objectiveType << std::endl;
    exit(-1);
  }
}

objectiveFn getObjective(GAParams * p) {
  if (p->objective.compare("acceptorPeak") == 0) {
    return objAcceptorPeak;
  }
  else if (p->objective.compare("acceptorAvg") == 0) {
    return objAcceptorAvg;
  }
  else if (p->objective.compare("acceptorAvgAfterPeak") == 0) {
    return objAcceptorAvgAfterPeak;
  }
  else if (p->objective.compare("acceptorFinal") == 0) {
    return objAcceptorFinal;
  }
  else {
    std::cerr << "[" << __FUNCTION__ << "]: unrecognized objective" << p->objective << std::endl;
    exit(-1);
  }
}