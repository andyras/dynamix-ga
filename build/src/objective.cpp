#include "objective.hpp"

// #define DEBUG

float singleObjective(GAGenome &c) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;
  // get MPI rank and size
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

#ifdef DEBUG
  std::cout << "rank: " << rank << " of " << size << "\n";
#endif

  // output value
  float output = 1.0;

  // Struct of parameters //////////////////////////////////////////////////////

  GAParams gp;
  assignGAParams("ins/ga.in", &gp);
  Params * p = &(gp.p);
  // TODO XXX FUGLINESS HARD-CODING
  std::string jobPrefix = "./";
  p->inputFile = jobPrefix + "ins/parameters.in";
  p->cEnergiesInput = jobPrefix + "ins/c_energies.in";
  p->bEnergiesInput = jobPrefix + "ins/b_energies.in";
  p->VNoBridgeInput = jobPrefix + "ins/Vnobridge.in";
  p->VBridgeInput = jobPrefix + "ins/Vbridge.in";

  // assign parameters from input file /////////////////////////////////////////

  assignParams(p->inputFile.c_str(), p);

  // assign GA parameters
  // this function is independent of the objective you are using. It determines
  // what the relevant parameters are for the optimization.
  if (gp.variables.compare("g1g2g1_c") == 0) {
    init_gammas(c, p);
  }
  else if (gp.variables.compare("wavepacket") == 0) {
    init_wavepacket(c, p);
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "variable set" <<
      gp.variables << "not recognized." << std::endl;
    exit(-1);
  }
  print1DGenes(c);

  initialize(p);

  // set number of processors for OpenMP ///////////////////////////////////////

  omp_set_num_threads(p->nproc);
  mkl_set_num_threads(p->nproc);

  // Make plot files ///////////////////////////////////////////////////////////

  makePlots(p);

  // propagate /////////////////////////////////////////////////////////////////

  propagate(p);

  // calculate value of objective function /////////////////////////////////////
#ifdef DEBUG
  std::cout << "Calculating value of objective function..." << std::endl;
#endif
  // set 'output' according to an objective function ///////////////////////////
  if (gp.objective.compare("acceptorPeak") == 0) {
    output = objAcceptorPeak(p);
  }
  else if (gp.objective.compare("acceptorAvg") == 0) {
    output = objAcceptorAvg(p);
  }
  else if (gp.objective.compare("acceptorAvgAfterPeak") == 0) {
    output = objAcceptorAvgAfterPeak(p);
  }
  else if (gp.objective.compare("acceptorFinal") == 0) {
    output = objAcceptorFinal(p);
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "objective function" <<
      gp.objective << "not recognized." << std::endl;
    exit(-1);
  }

  int pid = getpid();
  std::cout << "[" << pid << ":" << rank << "] " << "Objective: " << output << std::endl;


  return output;
}

float doubleObjective(GAGenome &c) {
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

  GAParams gp;
  assignGAParams("ins/ga.in", &gp);
  Params * p = &(gp.p);

  // TODO XXX FUGLINESS HARD-CODING
  std::string jobPrefix = "./";
  p->inputFile = jobPrefix + "ins/parameters.in";
  p->cEnergiesInput = jobPrefix + "ins/c_energies.in";
  p->bEnergiesInput = jobPrefix + "ins/b_energies.in";
  p->VNoBridgeInput = jobPrefix + "ins/Vnobridge.in";
  p->VBridgeInput = jobPrefix + "ins/Vbridge.in";

  // coherent propagation //////////////////////////////////////////////////////

  // assign parameters from input file /////////////////////////////////////////

  assignParams(p->inputFile.c_str(), p);
  p->coherent = 1;

  // assign GA parameters
  // this function is independent of the objective you are using. It determines
  // what the relevant parameters are for the optimization.
  if (gp.variables.compare("g1g2g1_c") == 0) {
    init_gammas(c, p);
  }
  else if (gp.variables.compare("wavepacket") == 0) {
    init_wavepacket(c, p);
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "variable set" <<
      gp.variables << "not recognized." << std::endl;
    exit(-1);
  }

  initialize(p);

  print1DGenes(c);

  // Make plot files ///////////////////////////////////////////////////////////

  makePlots(p);

  // set number of processors for OpenMP ///////////////////////////////////////

  omp_set_num_threads(p->nproc);
  mkl_set_num_threads(p->nproc);

  // propagate /////////////////////////////////////////////////////////////////

  propagate(p);

  // calculate value of objective function /////////////////////////////////////
#ifdef DEBUG
  std::cout << "Calculating value of objective function..." << std::endl;
#endif
  // set 'output' according to an objective function ///////////////////////////
  if (gp.objective.compare("acceptorPeak") == 0) {
    output1 = objAcceptorPeak(p);
  }
  else if (gp.objective.compare("acceptorAvg") == 0) {
    output1 = objAcceptorAvg(p);
  }
  else if (gp.objective.compare("acceptorAvgAfterPeak") == 0) {
    output1 = objAcceptorAvgAfterPeak(p);
  }
  else if (gp.objective.compare("acceptorFinal") == 0) {
    output1 = objAcceptorFinal(p);
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "objective function" <<
      gp.objective << "not recognized." << std::endl;
    exit(-1);
  }

  std::cout << "[" << pid << ":" << rank << "] " << "Objective 1: " << output1 << std::endl;

  // incoherent propagation ////////////////////////////////////////////////////
  assignParams(p->inputFile.c_str(), p);
  p->coherent = 0;
  if (gp.variables.compare("g1g2g1_c") == 0) {
    init_gammas(c, p);
  }
  else if (gp.variables.compare("wavepacket") == 0) {
    init_wavepacket(c, p);
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "variable set" <<
      gp.variables << "not recognized." << std::endl;
    exit(-1);
  }
  initialize(p);
  print1DGenes(c);
  propagate(p);
#ifdef DEBUG
  std::cout << "Calculating value of objective function again..." << std::endl;
#endif
  if (gp.objective.compare("acceptorPeak") == 0) {
    output2 = objAcceptorPeak(p);
  }
  else if (gp.objective.compare("acceptorAvg") == 0) {
    output2 = objAcceptorAvg(p);
  }
  else if (gp.objective.compare("acceptorAvgAfterPeak") == 0) {
    output2 = objAcceptorAvgAfterPeak(p);
  }
  else if (gp.objective.compare("acceptorFinal") == 0) {
    output2 = objAcceptorFinal(p);
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "objective function" <<
      gp.objective << "not recognized." << std::endl;
    exit(-1);
  }
  std::cout << "[" << pid << ":" << rank << "] " << "Objective 2: " << output2 << std::endl;

  double f = fabs(output1 - output2);
  std::cout << "[" << pid << ":" << rank << "] " << "Combined objective: " << f << std::endl;

  return f;
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

  f = 1.0 - findArrayMaximum(&(tcprobs[0]), tcprobs.size())/summ;
#ifdef DEBUG
  std::cout << "Value of f is 1.0 - " << findArrayMaximum(&(tcprobs[0]), p->numOutputSteps) << "/" << summ << " = " << f << std::endl <<
    "sigma: " << p->bulkGaussSigma << " mu: " << p->bulkGaussMu << " f: " << f << std::endl;
#endif

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

  f = 1.0 - tcprob/summ/(p->numOutputSteps + 1);

#ifdef DEBUG
  std::cout << "Value of f is 1.0 - " << tcprob << "/" << summ << "/" << (p->numOutputSteps + 1) << " = " << f << std::endl <<
    "sigma: " << p->bulkGaussSigma << " mu: " << p->bulkGaussMu << " f: " << f << std::endl;
#endif

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

  f = 1.0 - tcprob/summ/(p->numOutputSteps + 1 - maxInd);

#ifdef DEBUG
  std::cout << "Value of f is 1.0 - " << tcprob << "/" << summ << "/" << (p->numOutputSteps + 1 - maxInd) << " = " << f << std::endl <<
    "sigma: " << p->bulkGaussSigma << " mu: " << p->bulkGaussMu << " f: " << f << std::endl;
#endif

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

  f = 1.0 - finalPop;

#ifdef DEBUG
  std::cout << "sigma: " << p->bulkGaussSigma << " mu: " << p->bulkGaussMu << " f: " << f << std::endl;
#endif

  return f;
}
