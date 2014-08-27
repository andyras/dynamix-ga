#include <boost/mpi.hpp>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <dynamix.hpp>
#include <ga-mpi/ga.h>
#include <ga-mpi/std_stream.h>
#include <params.hpp>

#include "initializeGenome.hpp"
#include "mutator.hpp"
#include "output.hpp"
#include "parser.hpp"

// #define DEBUG

int mpi_tasks, mpi_rank;

int main(int argc, char **argv)
{
  // MPI init
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;

  int pid = getpid();

  // See if we've been given a seed to use (for testing purposes).  When you
  // specify a random seed, the evolution will be exactly the same each time
  // you use that seed number
  unsigned int seed = 0;
  for(int i=1 ; i<argc ; i++)
    if(strcmp(argv[i++],"seed") == 0)
      seed = atoi(argv[i]);

  // Declare variables for the GA parameters and set them to some default values.
  int ngen = 200; // Generations

  GAParams gp;
  dynamixGAParams dgp;

  if (world.rank() == 0) {
    assignGAParams("./ins/ga.in", &dgp);
    assignParams("./ins/parameters.in", &(dgp.gp.p));
    initialize(&(dgp.gp.p));
    for (int ii = 1; ii < world.size(); ii++) {
      world.send(ii, 0, dgp);
    }
  }
  else {
    world.recv(0, 0, dgp);
  }

  void * userData = &dgp;
#ifdef DEBUG
  std::cout << "[" << pid << ":" << mpi_rank <<  "] userData is " << userData << std::endl;
#endif

#ifdef DEBUG
  std::cout << "[" << pid << ":" << mpi_rank <<  "] Using " << gp.objectiveType << " objective function." << std::endl;
#endif

  unsigned int genomeLength = 0;
  if (dgp.gp.initializer.compare("g1g2g1_c") == 0) {
    genomeLength = 3;
    dgp.gp.lb.resize(genomeLength);
    dgp.gp.ub.resize(genomeLength);
    for (unsigned int ii = 0; ii < 3; ii++) {
      dgp.gp.lb[ii] = 0.001;
      dgp.gp.ub[ii] = 0.05;
    }
  }
  else if (dgp.gp.initializer.compare("wavepacket") == 0) {
    genomeLength = 2;
    dgp.gp.lb.resize(genomeLength);
    dgp.gp.ub.resize(genomeLength);
    for (unsigned int ii = 0; ii < 2; ii++) {
      dgp.gp.lb[ii] = 0.0001; // this needs to be ~> the interlevel spacing
      dgp.gp.ub[ii] = 0.01; // it would be nice to have these initialized by reading parameters.in
    }
  }
  else if (dgp.gp.initializer.compare("wavepacketGammas") == 0) {
    genomeLength = 5;
    dgp.gp.lb.resize(genomeLength);
    dgp.gp.ub.resize(genomeLength);
    for (unsigned int ii = 0; ii < 3; ii++) {
      dgp.gp.lb[ii] = 0.001;
      dgp.gp.ub[ii] = 0.05;
    }
    for (unsigned int ii = 3; ii < (3+2); ii++) {
      dgp.gp.lb[ii] = 0.0001;
      dgp.gp.ub[ii] = 0.01;
    }
  }
  else if (dgp.gp.initializer.compare("torsion") == 0) {
    genomeLength = 3;
    dgp.gp.lb.resize(genomeLength);
    dgp.gp.ub.resize(genomeLength);
    // torsionSin2V0
    dgp.gp.lb[0] = 0.00001;
    dgp.gp.ub[0] = 0.001;
    // torsionSin2V1
    dgp.gp.lb[1] = 0.00001;
    dgp.gp.ub[1] = 0.001;
    // Eb
    dgp.gp.lb[2] = dgp.gp.p.kBandEdge;
    dgp.gp.ub[2] = dgp.gp.p.kBandTop;
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "variable set" <<
      dgp.gp.initializer << "not recognized." << std::endl;
    exit(-1);
  }

  GA1DArrayGenome<double> genome(genomeLength, dgp.getObjectiveType(), userData);
  genome.initializer(dgp.getInitializer());

  genome.mutator(getMutator(&(dgp.gp)));

  // initialize structures for best genomes/scores
  double initVal = 0.0;
  if (dgp.gp.minmax.compare("min") == 0) {
    initVal = INFINITY;
  }
  else if (dgp.gp.minmax.compare("max") == 0) {
    initVal = -INFINITY;
  }

  dgp.gp.bestScore = initVal;
  dgp.gp.bestGenome.resize(genomeLength);
  // Now create the GA using the genome and run it. We'll use sigma truncation
  // scaling so that we can handle negative objective scores.
  GASimpleGA ga(genome); // TODO change to steady-state
  ga.userData(userData);
  GALinearScaling scaling;
  if (dgp.gp.minmax.compare("min") == 0) {
    ga.minimize();
    ga.pConvergence(1.0 + dgp.gp.convergence);
  }
  else if (dgp.gp.minmax.compare("max") == 0) {
    ga.maximize();
    ga.pConvergence(1.0 - dgp.gp.convergence);
  }
  else {
    std::cerr << "ERROR: unrecognized minmax: << " << dgp.gp.minmax << std::endl;
    exit(0);
  }
  ga.populationSize(dgp.gp.popsize);
  ga.nGenerations(ngen);
  ga.pMutation(dgp.gp.pMut);
  ga.pCrossover(dgp.gp.pCross);
  ga.terminator(GAGeneticAlgorithm::TerminateUponConvergence);
  ga.scaling(scaling);
  if(mpi_rank == 0)
    ga.scoreFilename("evolution.txt");
  else
    ga.scoreFilename("/dev/null");
  ga.scoreFrequency(1);
  ga.flushFrequency(1);
  ga.selectScores(GAStatistics::AllScores);
  // Pass MPI data to the GA class
  ga.mpi_rank(mpi_rank);
  ga.mpi_tasks(mpi_tasks);
  ga.evolve(seed);

  std::vector<double> bestScores;
  std::vector<double> bestGenomes;

  if (mpi_rank == 0) {
    bestScores.resize(mpi_tasks, 0.0);
    bestGenomes.resize(mpi_tasks*genomeLength, 0.0);
  }

  // wait for everyone to finish ///////////////////////////////////////////////
  MPI_Barrier(MPI_COMM_WORLD);

  // gather best scores ////////////////////////////////////////////////////////
  MPI_Gather(&(dgp.gp.bestScore), 1, MPI_DOUBLE,
             &(bestScores[0]), 1, MPI_DOUBLE,
             0, MPI_COMM_WORLD);

  // gather best genomes ///////////////////////////////////////////////////////
  MPI_Gather(&(dgp.gp.bestGenome[0]), genomeLength, MPI_DOUBLE,
             &(bestGenomes[0]), genomeLength, MPI_DOUBLE,
             0, MPI_COMM_WORLD);

  // Dump the GA results to file ///////////////////////////////////////////////
  if(mpi_rank == 0)
  {
    // find rank of best score
    auto bestScore = std::min_element(bestScores.begin(), bestScores.end());
    if (dgp.gp.minmax.compare("max") == 0) {
      bestScore = std::max_element(bestScores.begin(), bestScores.end());
    }
    unsigned int bestIdx = std::distance(bestScores.begin(), bestScore);

    std::cout << "[" << pid << ":" << mpi_rank << "] Best score: " <<
      bestScores[bestIdx] << " Genome:";
    for (unsigned int ii = 0; ii < genomeLength; ii++) {
      std::cout << " " << bestGenomes[bestIdx*genomeLength + ii];
    }
    std::cout << std::endl;

    std::cout << "[" << pid << ":" << mpi_rank << "] Scores:";
    for (unsigned int ii = 0; ii < bestScores.size(); ii++) {
      std::cout << " " << bestScores[ii];
    }
    std::cout << std::endl;

    // create new genome and evaluate objective (in order to create outputs)
    GA1DArrayGenome<double> bg(genomeLength, dgp.getObjectiveType(), userData);
    // fill genome with best parameters
    for (int ii = 0; ii < bg.length(); ii++) {
      bg.gene(ii, bestGenomes[bestIdx*genomeLength + ii]);
    }
    // call objective function
    (* dgp.getObjectiveType())(bg);
  }

  MPI_Finalize();

  return 0;
}
