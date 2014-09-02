#include "objective.hpp"

// #define DEBUG

float singleObjective(GAGenome &c) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;
  dynamixGAParams * dgp = (dynamixGAParams *) genome.userData();
  GAParams * gp = &(dgp->gp);

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
  initParamsToChange(c, &p);

  initialize(&p);
  makePlots(&p);
  propagate(&p);

  // calculate value of objective function /////////////////////////////////////
#ifdef DEBUG
  std::cout << "[" << pid << ":" << rank << "] " <<
    "Calculating value of objective function..." << std::endl;
#endif
  output = 0.0;
  std::vector< std::tuple<objectiveFn, double> > obs = dgp->ot.getObjWeights();
  objectiveFn obj; // objective function to call
  double w; // weight
  // call each objective fn and add the weighted value to the output ///////////
  for (unsigned int ii = 0; ii < obs.size(); ii++) {
    obj = std::get<0>(obs[ii]);
    w = std::get<1>(obs[ii]);
#ifdef DEBUG
    std::cout << "objective value is " << (*obj)(&p) << " and weight is " << w << std::endl;
#endif
    output += (*obj)(&p)*w; // dereference the function pointer and call it
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
  dynamixGAParams * dgp = (dynamixGAParams *) genome.userData();
  GAParams * gp = &(dgp->gp);

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
  initParamsToChange(c, &p);

  // coherent propagation //////////////////////////////////////////////////////

  p.coherent = 1;
  p.outputDir = "./coh/";

  initialize(&p);
  makePlots(&p);
  propagate(&p);

  // calculate value of objective function /////////////////////////////////////
#ifdef DEBUG
  std::cout << "[" << pid << ":" << rank << "] " <<
    "Calculating value of objective function..." << std::endl;
#endif
  output1 = 0.0;
  std::vector< std::tuple<objectiveFn, double> > obs = dgp->ot.getObjWeights();
  objectiveFn obj; // objective function to call
  double w; // weight
  // call each objective fn and add the weighted value to the output ///////////
  for (unsigned int ii = 0; ii < obs.size(); ii++) {
    obj = std::get<0>(obs[ii]);
    w = std::get<1>(obs[ii]);
#ifdef DEBUG
    std::cout << "objective value is " << (*obj)(&p) << " and weight is " << w << std::endl;
#endif
    output1 += (*obj)(&p)*w; // dereference the function pointer and call it
  }

  std::cout << "[" << pid << ":" << rank << "] " << "Objective 1: " << output1 << " Genome:";
  for (unsigned int ii = 0; ii < (unsigned int) genome.length(); ii++) {
    std::cout << " " << genome.gene(ii);
  }
  std::cout << std::endl;

  // incoherent propagation ////////////////////////////////////////////////////
  p.coherent = 0;
  p.outputDir = "./inc/";

  initialize(&p);
  makePlots(&p);
  propagate(&p);

  // calculate value of objective function /////////////////////////////////////
#ifdef DEBUG
  std::cout << "[" << pid << ":" << rank << "] " <<
    "Calculating value of objective function again..." << std::endl;
#endif
  output2 = 0.0;
  // call each objective fn and add the weighted value to the output ///////////
  for (unsigned int ii = 0; ii < obs.size(); ii++) {
    obj = std::get<0>(obs[ii]);
    w = std::get<1>(obs[ii]);
#ifdef DEBUG
    std::cout << "objective value is " << (*obj)(&p) << " and weight is " << w << std::endl;
#endif
    output2 += (*obj)(&p)*w; // dereference the function pointer and call it
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